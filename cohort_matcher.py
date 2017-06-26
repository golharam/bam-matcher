#!/bin/env python
'''
cohort_matcher - Compare genotypes of two cohorts of samples
'''
import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys
import logging
import pysam
import vcf
from fisher import pvalue

logger = logging.getLogger(__name__)

def checkConfig(config):
    '''
    checkConfig makes sure all the parameters specified are valid
    '''
    if os.path.isdir(config.cache_dir) is False:
        logger.error("Cache dir (%s) is not a directory", config.cache_dir)
        return False
    if os.path.isdir(config.scratch_dir) is False:
        logger.error("Scratch dir (%s) is not a directory", config.scratch_dir)
        return False
    if os.path.isdir(config.output_dir) is False:
        logger.error("Output directory (%s) does not exist", config.output_dir)
        return False

    if os.path.exists(config.vcf) is False:
        logger.error("VCF file (%s) is not accessible", config.vcf)
        return False

    if os.path.exists(config.reference) is False:
        logger.error("Reference FASTA file (%s) is not accessible", config.reference)
        return False
    reference_index = config.reference + ".fai"
    if os.path.exists(reference_index) is False:
        logger.error("Reference FASTA file (%s) is not indexed", reference_index)
        return False

    # If reference2, vcf2, or chromosome-map is specified, make sure all are specified and exist
    if config.reference2 is not None or config.vcf2 is not None or config.chromosome_map is not None:

        if config.reference2 is None:
            logger.error("Reference2 must be specified if vcf2 or chromosome-map is specified")
            return False
        elif config.vcf2 is None:
            logger.error("vcf2 must be specified if reference2 or chromosome-map is specified")
            return False
        elif config.chromosome_map is None:
            logger.error("Chromosome map must be specified is reference2 or vcf2 is specified")
            return False

        if os.path.exists(config.vcf2) is False:
            logger.error("VCF2 (%s) is not accessible", config.vcf2)
            return False

        if os.path.exists(config.reference2) is False:
            logger.error("Reference2 FASTA file (%s) is not accessible", config.reference2)
            return False
        reference2_index = config.reference2 + ".fai"
        if os.path.exists(reference2_index) is False:
            logger.error("Reference 2 FASTA file (%s) is not indexed", reference2_index)
            return False

        if os.path.exists(config.chromosome_map) is False:
            logger.error("Chromosome map (%s) could not be read", config.chromosome_map)
            return False

    if config.caller == "freebayes" and config.freebayes_path is None:
        logger.error("Freebayes-path must be set when caller is Freebayes")
        return False

    if os.path.exists(config.aws) is False:
        logger.error("Unable to locate AWS CLI: %s", config.aws)
        return False

    if os.path.exists(config.freebayes_path) is False:
        logger.error("Unable to locate caller: %s", config.freebayes_path)
        return False

    if os.path.exists(config.samtools) is False:
        logger.error("Unable to locate samtools: %s", config.samtools)
        return False

    return True

def checkReference(sample, localBamFile, reference, vcf):
    ''' Make sure the reference used BAM file is the same as
        the vcf file
    '''
    bam_chroms = get_chrom_names_from_BAM(localBamFile)
    logger.debug("BAM chromosomes: %s", bam_chroms)

    REF_CHROMS = get_chrom_names_from_REF(reference)
    logger.debug("REF chromosomes: %s", REF_CHROMS)

    if set(REF_CHROMS).issubset(set(bam_chroms)) == False:
        bamREF_dff = set(REF_CHROMS).difference(set(bam_chroms))
        logger.error("Sample BAM %s contains chromosomes not in reference %s: %s",
                     sample, reference, bamREF_dff)
        return False

    ''' Make sure the VCF file and the BAM file have matching chromosome names '''
    VCF_chroms = get_chrom_names_from_VCF(vcf)
    logger.debug("VCF chromosomes: %s", VCF_chroms)

    if set(VCF_chroms).issubset(set(bam_chroms)) == False:
        bamVCF_diff = set(VCF_chroms).difference(set(bam_chroms))
        logger.error("Sample VCF %s contains chromosomes not in BAM %s: %s",
                     sample, vcf, bamVCF_diff)
        return False
    return True

def compareGenotypes(var_list, var_list2, intersection, alternate_chroms, def_to_alt):
    ct_common = 0
    comm_hom_ct = 0
    comm_het_ct = 0
    ct_diff = 0
    diff_hom_ct = 0
    diff_het_ct = 0
    diff_1sub2_ct = 0
    diff_hom_het_ct = 0
    diff_2sub1_ct = 0
    diff_het_hom_ct = 0
    for pos_ in intersection:
        gt1 = var_list[pos_]['GT']
        if alternate_chroms is not None:
            gt2 = def_to_alt[pos_]
        else:
            gt2 = var_list2[pos_]['GT']
        #gt1 = bam1_gt[pos_]
        #gt2 = bam2_gt[pos_]
        # if genotypes are the same
        if is_same_gt(gt1, gt2):
            ct_common += 1
            if is_hom(gt1):
                comm_hom_ct += 1
            else:
                comm_het_ct += 1
        else:
            ct_diff += 1
            # both are hom and different
            if is_hom(gt1) and is_hom(gt2):
                diff_hom_ct += 1
            # both are het and different
            elif is_hom(gt1) == False and is_hom(gt2) == False:
                diff_het_ct += 1
            # one is hom, one is het, test for subset
            elif is_hom(gt1):
                if is_subset(gt1, gt2):
                    diff_1sub2_ct += 1
                else:
                    diff_hom_het_ct += 1
            elif is_hom(gt2):
                if is_subset(gt2, gt1):
                    diff_2sub1_ct += 1
                else:
                    diff_het_hom_ct += 1
            else:
                print "WTF?"
                print gt1, gt2

    total_compared = ct_common + ct_diff
    frac_common_plus = 0
    frac_common = 0
    if total_compared > 0:
        frac_common = float(ct_common)/total_compared
        frac_common_plus = float(ct_common + max(diff_2sub1_ct,
                                                 diff_1sub2_ct))/total_compared

    # test of allele-specific genotype subsets
    allele_subset = ""
    sub_sum = diff_1sub2_ct + diff_2sub1_ct
    # don't bother if fewer than 10
    if sub_sum > 10:
        pv_set = pvalue(diff_1sub2_ct, diff_2sub1_ct, ct_diff/2, ct_diff/2)
        pv_ = min(pv_set.left_tail, pv_set.right_tail)
        if pv_ < 0.05:
            if diff_1sub2_ct < diff_2sub1_ct:
                allele_subset = "2sub1"
                frac_common_plus = float(ct_common + diff_2sub1_ct) / total_compared
            else:
                allele_subset = "1sub2"
                frac_common_plus = float(ct_common + diff_1sub2_ct) / total_compared

    results = {}
    results['total_compared'] = total_compared
    results['ct_common'] = ct_common
    results['frac_common'] = frac_common
    results['frac_common_plus'] = frac_common_plus
    results['comm_hom_ct'] = comm_hom_ct
    results['comm_het_ct'] = comm_het_ct
    results['ct_diff'] = ct_diff
    results['diff_hom_ct'] = diff_hom_ct
    results['diff_het_ct'] = diff_het_ct
    results['diff_hom_het_ct'] = diff_hom_het_ct
    results['diff_het_hom_ct'] = diff_het_hom_ct
    results['diff_1sub2_ct'] = diff_1sub2_ct
    results['diff_2sub1_ct'] = diff_2sub1_ct
    results['allele_subset'] = allele_subset
    results['judgement'], results['short_judgement'] = makeJudgement(total_compared, frac_common,
                                                                     frac_common_plus, allele_subset)
    return results

def compareSamples(sampleSet1, sampleSet2, config):
    ''' Compare all the samples against each other '''
    if config.chromosome_map:
        default_chroms, alternate_chroms, def_to_alt, alt_to_def = get_chrom_data_from_map(config.chromosome_map)
    else:
        default_chroms, alternate_chroms, def_to_alt, alt_to_def = None, None, None, None

    frac_common_matrix = {}
    total_compared_matrix = {}
    for sample1 in sampleSet1:
        for sample2 in sampleSet2:
            logger.info("Comparing {} - {}".format(sample1["name"], sample2["name"]))
            ''' Get a list of variants that pass in sample 1 '''
            tsv1 = os.path.join(config.cache_dir, sample1["name"] + ".tsv")
            var_list = get_tsv_variants(tsv1, config.dp_threshold)
            ''' then parse second tsv file to get list of variants that passed in both samples '''
            tsv2 = os.path.join(config.cache_dir, sample2["name"] + ".tsv")
            var_list2 = get_tsv_variants(tsv2, config.dp_threshold)

            intersection = getIntersectingVariants(var_list, var_list2, default_chroms, alternate_chroms)

            ''' compare the genotypes '''
            results = compareGenotypes(var_list, var_list2, intersection, alternate_chroms, def_to_alt)
            writeSampleComparisonReport(sample1["name"], sample2["name"], config, results)
            
            ''' save to grand matrix '''
            if sample1["name"] not in frac_common_matrix:
                frac_common_matrix[sample1["name"]] = {}
                total_compared_matrix[sample1["name"]] = {}
            frac_common_matrix[sample1["name"]][sample2["name"]] = results['frac_common']
            total_compared_matrix[sample1["name"]][sample2["name"]] = results['total_compared']
            
    writeSimilarityMatrix(config, sampleSet1, sampleSet2, frac_common_matrix, total_compared_matrix)
    
    
def downloadBAMFile(bamFile, config):
    localBamFile = os.path.join(config.scratch_dir, os.path.basename(bamFile))
    if os.path.exists(localBamFile):
        logger.debug("Using cached bam file: {}".format(localBamFile))
    else:
        if downloadFileFromAmazon(bamFile, config.scratch_dir, config) is None:
            logger.error("File does not exist or is inaccessible in Amazon.")
            return None

    ''' If the index is already downloaded, use it '''
    bamIndex1 = bamFile.rstrip(".bam") + ".bai"
    bamIndex2 = bamFile + ".bai"
    localBamIndex1 = localBamFile.rstrip(".bam") + ".bai"
    localBamIndex2 = localBamFile + ".bai"
    localBamIndex = ''

    if os.path.exists(localBamIndex1):
        localBamIndex = localBamIndex1
    elif os.path.exists(localBamIndex2):
        localBamIndex = localBamIndex2
    else:
        ''' Else, try to download index '''
        if downloadFileFromAmazon(bamIndex1, config.scratch_dir, config):
            localBamIndex = localBamIndex1
        else:
            ''' Try to download index 2 '''
            if downloadFileFromAmazon(bamIndex2, config.scratch_dir, config):
                localBamIndex = localBamIndex2
            else:
                logger.debug("Could not find matching bam index.  Generating...")
                cmd = [config.samtools, 'index', localBamFile]
                p = subprocess.Popen(cmd)
                p.wait()
                if p.returncode != 0:
                    logger.error("Unable to generated BAM index")
                    exit(1)
                localBamIndex = localBamIndex2
    return localBamFile, localBamIndex

def downloadFileFromAmazon(srcFile, destDirectory, config):
    if isFileInAmazon(srcFile, config) == False:
        return None

    cmd = [config.aws, "s3", "cp", srcFile, destDirectory]
    logger.info("Downloading file: {}".format(srcFile))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()

    if p.returncode != 0:
        return None

    localFile = os.path.join(destDirectory, os.path.basename(srcFile))
    if os.access(localFile, os.R_OK) == False:
        return None

    return localFile

def get_chrom_names_from_BAM(bam_file):
    ''' Get list of chromosome names from the BAM file '''
    chrom_list = []
    inbam = pysam.AlignmentFile(bam_file, "rb")
    header_sq = inbam.header["SQ"]
    for sq_ in header_sq:
        chrom_list.append(sq_["SN"])
    return chrom_list

def get_chrom_names_from_REF(ref_fasta):
    ''' Get list of chromosome names from the reference file '''
    ref_idx = ref_fasta + ".fai"
    chrom_list = []
    fin = open(ref_idx, "r")
    for line in fin:
        if line.strip() == "":
            continue
        chrom_list.append(line.strip().split()[0])
    return chrom_list

def get_chrom_data_from_map(chrom_map_file):
    chrom_ct = 0
    default_chroms = []
    alternate_chroms = []
    def_to_alt = {}
    alt_to_def = {}

    with open(chrom_map_file, 'r') as fin:
        header = fin.readline()
        for line in fin:
            if line.strip() == "":
                continue
            chrom_ct += 1
            bits = line.strip().split()
            default_chroms.append(bits[0])
            alternate_chroms.append(bits[1])
            def_to_alt[bits[0]] = bits[1]
            alt_to_def[bits[1]] = bits[0]
    return default_chroms, alternate_chroms, def_to_alt, alt_to_def

def get_chrom_names_from_VCF(vcf_file):
    ''' Get list of chromosome names from the VCF file '''
    chrom_list = []
    with open(vcf_file, "r") as fin:
        vcf_reader = vcf.Reader(fin)
        for vcfRecord in vcf_reader:
            if chrom_list.count(vcfRecord.CHROM) == 0:
                chrom_list.append(vcfRecord.CHROM)
    return chrom_list

def getIntersectingVariants(var_list, var_list2, alternate_chroms, alt_to_def):
    ''' Given two tsv list of variants, find the intersection '''
    intersection = []
    for var2 in var_list2:
        bits = var2.split("\t")
        if alternate_chroms is not None:
            var_ = alt_to_def[bits[0]] + "\t" + bits[1]
        else:
            var_ = var2
        
        if var_ in var_list:
            intersection.append(var_)

    return intersection

def get_tsv_variants(tsvFile, dp_threshold):
    ''' Return a list of variants that pass threshold '''
    var_list = {}
    with open(tsvFile, "r") as fin:
        for line in fin:
            if line.startswith("CHROM\t"):
                continue
            bits = line.strip("\n").split("\t")
            if bits[5] == "NA":
                continue
            elif int(bits[5]) < dp_threshold:
                continue
            else:
                var_list["\t".join(bits[:2])] = {'REF': bits[2], 'ALT': bits[3],
                                                 'DP': int(bits[5]), 'GT': bits[7]}
    return var_list

def genotypeSample(sample, bamFile, reference, vcf, intervalsFile, config):
    '''
    This function calls the actual genotype to generate a vcf, then 
    converts the vcf to tsv file.
    '''
    logger.info("Genotyping {}".format(sample))

    deleteBam = False
    outputVcf = os.path.join(config.cache_dir, sample + ".vcf")
    if os.path.exists(outputVcf) is False:
        # Download the BAM file and index if they are not local
        if bamFile.startswith("s3://"):
            localBamFile, localBamIndex = downloadBAMFile(bamFile, config)
            deleteBam = True
        else:
            localBamFile = os.path.abspath(bamFile)

        # Make sure BAM and reference have matching chromosomes
        logger.debug("Checking %s against reference", sample)
        if checkReference(sample, localBamFile, reference, vcf) is False:
            exit(1)

        if config.caller == 'freebayes':
            logger.debug("{}: Calling freebayes".format(sample))
            cmd = [config.freebayes_path, "--fasta-reference", reference, "--targets",
                   intervalsFile, "--no-indels", "--min-coverage", str(config.dp_threshold),
                   "--report-all-haplotype-alleles", "--report-monomorphic", "--vcf",
                   outputVcf, localBamFile]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p.wait()
            if p.returncode != 0:
                logger.error("Error executing %s.\nStdout: %s\nStderr: %s,", ' '.join(cmd), out, err)
                os.remove(outputVcf)
                exit(1)
            if os.path.exists(outputVcf) == False:
                logger.error("Output VCF file, {}, could not be found in cache.".format(outputVcf))
                exit(1)
        else:
            logger.error("other callers not yet supported")
            exit(1)
            
    # Convert the vcf to tsv
    logger.debug("{}: Convert VCF to TSV".format(sample))
    out_tsv = os.path.join(config.cache_dir, sample + ".tsv")
    if os.path.exists(out_tsv) == False:
        VCFtoTSV(outputVcf, out_tsv, config.caller)

    # Delete the bam/bai files if they were downloaded
    if deleteBam == True:
        os.remove(localBamFile)
        os.remove(localBamIndex)

    logger.info("{}: Done".format(sample))
    return None

def genotypeSamples(sampleSet, reference, vcf, intervalsFile, config):
    pool = multiprocessing.Pool(config.max_jobs)
    for sampleIndex in range(len(sampleSet)):
        sample = sampleSet[sampleIndex]
        tsvFile = os.path.join(config.cache_dir, sample["name"] + ".tsv")
        if os.path.exists(tsvFile) is False:
            if config.max_jobs == 1:
                genotypeSample(sample["name"], sample["bam"], reference, vcf, intervalsFile, config)
            else:
                pool.apply_async(genotypeSample, (sample["name"], sample["bam"], reference,
                                                  vcf, intervalsFile, config))
    pool.close()
    pool.join()

def is_hom(gt):
    gt_ = gt.split("/")
    if gt_[0] == gt_[1]:
        return True
    else:
        return False

def is_same_gt(gt1, gt2):
    if gt1 == gt2:
        return True
    else:
        gt1_ = sorted(gt1.split("/"))
        gt2_ = sorted(gt2.split("/"))
        if gt1_ == gt2_:
            return True
        else:
            return False

def is_subset(hom_gt, het_gt):
    gt_hom = hom_gt.split("/")[0]
    if gt_hom in het_gt:
        return True
    else:
        return False

def isFileInAmazon(srcFile, config):
    cmd = [config.aws, "s3", "ls", srcFile]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()
    lines = out.splitlines()
    for line in lines:
        fields = line.split()
        if fields[3] == os.path.basename(srcFile):
            return True
    return False

def main(argv):
    config = parseArguments(argv)
    logging.basicConfig(level=config.log_level)
    logger.info("cohort-matcher")
    logger.info(config)

    if checkConfig(config) is False:
        return 1

    sampleSet1 = readSamples(config.set1)
    sampleSet2 = readSamples(config.set2)
    if sampleSet1 is False or sampleSet2 is False:
        return 1

    intervalsFile = os.path.join(config.cache_dir, "set1.intervals")
    intervalsFile2 = os.path.join(config.cache_dir, "set2.intervals")
    logger.info("Generating intervals file (%s) from VCF (%s)", intervalsFile, config.vcf)
    vcfToIntervals(config.vcf, intervalsFile)
    if config.reference2 is None:
        logger.info("Copying intervals file %s -> %s", intervalsFile, intervalsFile2)
        shutil.copyfile(intervalsFile, intervalsFile2)
    else:
        logger.info("Generating intervals file (%s) from VCF (%s)", intervalsFile2, config.vcf2)
        vcfToIntervals(config.vcf2, intervalsFile2)

    genotypeSamples(sampleSet1, config.reference, config.vcf, intervalsFile, config)
    if config.reference2 is None:
        genotypeSamples(sampleSet2, config.reference, config.vcf, intervalsFile2, config)
    else:
        genotypeSamples(sampleSet2, config.reference2, config.vcf2, intervalsFile2, config)

    compareSamples(sampleSet1, sampleSet2, config)
    return 0

def makeJudgement(total_compared, frac_common, frac_common_plus, allele_subset):
    ''' Make judgement of sample similarity based on genotype comparison '''
    A_BIT_LOW = """the number of comparable genomic loci is a bit low.
Try using a different variants list (--VCF) file which have more appropriate
genomic positions for comparison."""
    
    if total_compared <= 20:
        judgement = "Inconclusive: Too few loci to compare"
    elif total_compared <= 100:
        # allow for 0.90 frac_common for low loci count
        if frac_common >= 0.9 or frac_common_plus >= 0.9:
            judgement = "LIKELY SAME SOURCE: %s" % A_BIT_LOW
            short_judgement = "LIKELY SAME"
            # if there is possible allele-specific genotype subset
            if allele_subset == "1sub2" or allele_subset == "2sub1":
                sub_ = allele_subset.split("sub")[0]
                over_ = allele_subset.split("sub")[1]
                judgement += """BAM%s genotype appears to be a subset of BAM%s. Possibly
                BAM%s is RNA-seq data or BAM%s is contaminated.""" % (sub_, over_,
                                                                      sub_, over_)
                short_judgement += ". (BAM%s is subset of BAM%s)" % (sub_, over_)
        elif frac_common <= 0.6:
            judgement = "LIKELY FROM DIFFERENT SOURCES: %s" % A_BIT_LOW
            short_judgement = "LIKELY DIFFERENT"
        else:
            judgement = "INCONCLUSIVE: %s" % A_BIT_LOW
            short_judgement = "INCONCLUSIVE"
    # 3. >100 sites compared
    else:
        if frac_common >= 0.95:
            judgement = "BAM FILES ARE FROM THE SAME SOURCE"
            short_judgement = "SAME"
        elif frac_common_plus >= 0.95:
            short_judgement = "SAME"
            judgement = "BAM FILES ARE VERY LIKELY FROM THE SAME SOURCE"
            if allele_subset == "1sub2" or allele_subset == "2sub1":
                sub_ = allele_subset.split("sub")[0]
                over_ = allele_subset.split("sub")[1]
                judgement += """, but with possible allele specific genotype.\nBAM%s
                genotype appears to be a subset of BAM%s. Possibly BAM%s is RNA-seq
                data or BAM%s is contaminated. """ % (sub_, over_, sub_, over_)
                short_judgement += ". (BAM%s is subset of BAM%s)" % (sub_, over_)
        elif frac_common <= 0.6:
            judgement = "BAM FILES ARE FROM DIFFERENT SOURCES"
            short_judgement = "DIFFERENT"
        elif frac_common >= 0.8:
            judgement = """LIKELY FROM THE SAME SOURCE. However, the fraction of sites
            with common genotype is lower than expected. This can happen with samples
            with low coverage."""
            short_judgement = "LIKELY SAME"
        else:
            judgement = "LIKELY FROM DIFFERENT SOURCES"
            short_judgement = "LIKELY DIFFERENT"

    return judgement, short_judgement

def parseArguments(argv):
    parser = argparse.ArgumentParser(description="Compare two sets cohorts of bam files \
        to see if they are from the same samples, using frequently occuring SNPs \
        reported in the 1000Genome database")
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser_grp1 = parser.add_argument_group("Required")
    parser_grp1.add_argument("--set1", "-S1", required=True, help="First set of samples")
    parser_grp1.add_argument("--set2", "-S2", required=True, help="Second set of samples")

    parser_grp2 = parser.add_argument_group("Directories")
    parser_grp2.add_argument("--cache-dir", "-CD", required=False, default="./cache",
                             help="Specify directory for cached data. (Default: ./cache)")
    parser_grp2.add_argument("--scratch-dir", "-SD", required=False, default="/scratch",
                             help="Specify scratch directory. (Default: /scratch)")

    parser_grp3 = parser.add_argument_group("Genotyper")
    parser_grp3.add_argument("--caller", "-CL", required=False, default="freebayes",
                             choices=('gatk', 'freebayes', 'varscan'),
                             help="Specify which caller to use (default = 'freebayes')")
    parser_grp3.add_argument("--dp-threshold", "-DP", required=False, default=15, type=int,
                             help="Minimum required depth for comparing variants")
    parser_grp3.add_argument("--number_of_snps", "-N", required=False, type=int,
                             help="Number of SNPs to compare.")
    parser_grp3.add_argument("--max-jobs", required=False, default=1, type=int,
                             help="Maximum number of parallel genotyping jobs (Default: 1)")

    parser_grp4 = parser.add_argument_group("Reference")
    parser_grp4.add_argument("--vcf", "-V", required=True,
                             help="VCF file containing SNPs to check (default can be specified \
                                  in config file instead)")
    parser_grp4.add_argument("--vcf2", "-V2", default=None,
                             help="VCF file containing SNPs to check (default can be specified \
                                  in config file instead) when using reference2")
    parser_grp4.add_argument("--reference", "-R", required=True,
                             help="Reference FASTA File (indexed with samtools faidx) for set1 \
                             or both sets")
    parser_grp4.add_argument("--reference2", "-R2", default=None,
                             help="Reference FASTA File (indexed with samtools faidx) for set2")
    parser_grp4.add_argument("--chromosome-map", "-CM", default=None,
                             help="Chromosome mapping, if two reference are used")

    parser_grp5 = parser.add_argument_group("Freebayes Settings")
    parser_grp5.add_argument("--freebayes-path", required=False, default="freebayes",
                             help="Path to freebayes binary (if not in PATH")

    parser_grp7 = parser.add_argument_group("Paths")
    parser_grp7.add_argument("--aws", required=False, default="/usr/bin/aws",
                             help="Specify path to aws cli")
    parser_grp7.add_argument("--samtools", required=False, default="samtools",
                             help="Specify path to samtools")

    parser_grp8 = parser.add_argument_group("Output")
    parser_grp8.add_argument("--output-dir", "-O", required=False,
                             default="./cohort-matcher-output",
                             help="Specify output directory for sample comparisons \
                             (Default: ./cohort-matcher-output)")
    parser_grp8.add_argument("--short-output", "-so", required=False, default=False,
                             action="store_true", help="Short output format (Default: False")
    parser_grp8.add_argument("--report-file", required=False, default="cohort-matcher-results.txt",
                             help="Specify name of similarity matrix file \
                             (Default: ./cohort-master-results.txt)")
    '''
    # Freebayes options
    parser_grp3.add_argument("--fastfreebayes",  "-FF", required=False, default=False,
                             action="store_true", help="Use --targets option for Freebayes.")
    # GATK options
    parser_grp3.add_argument("--gatk-mem-gb" ,   "-GM", required=False,
                             type=int, help="Specify Java heap size for GATK (GB, int)")
    parser_grp3.add_argument("--gatk-nt" ,   "-GT", required=False,
                             type=int, help="Specify number of threads for GATK UnifiedGenotyper \
                             (-nt option)")
    # VarScan Options
    parser_grp3.add_argument("--varscan-mem-gb", "-VM", required=False,
                             type=int, help="Specify Java heap size for VarScan2 (GB, int)")
    '''
    args = parser.parse_args(argv)
    return args

def readSamples(sampleSheetFile):
    '''
    readSamples reads in a sampleSheetFile consisting of two columns:
    name and bamfile
    '''
    if os.path.isfile(sampleSheetFile) is False:
        logger.error("%s does not exist", sampleSheetFile)
        return False
    logger.info("Reading %s", sampleSheetFile)
    samples = []
    with open(sampleSheetFile, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if len(line) == 0 or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) != 2:
                logger.error("Expect 2 fields (sampleName, bamFile) but encountered %d",
                             len(fields))
                return False

            sample = {"name": fields[0],
                      "bam": fields[1]}
            samples.append(sample)
    logger.info("Read %d samples.", len(samples))
    return samples

def vcfToIntervals(vcfFile, bedFile, window=0, format="freebayes", cmap=None):
    '''
    Convert a vcf file to a 3-column interval/bed file
    '''
    if os.path.exists(bedFile) is True:
        logger.info("%s already exists.", bedFile)
        return

    vcf_read = vcf.Reader(open(vcfFile, "r"))
    fout = open(bedFile, "w")
    for var in vcf_read:
        # intervals format
        chrom_ = var.CHROM
        if cmap != None:
            if var.CHROM not in cmap:
                continue
            else:
                chrom_ = cmap[var.CHROM]
        if format == "gatk" or format == "varscan":
            start_pos = var.POS - window
            end_pos = start_pos + window
            fout.write("%s:%d-%d\n" % (chrom_, start_pos, end_pos))
        # BED format
        elif format == "bed" or format == "freebayes":
            start_pos = var.POS - window - 1
            end_pos = start_pos + window + 1
            fout.write("%s\t%d\t%d\n" % (chrom_, start_pos, end_pos))
    fout.close()
    return

def VCFtoTSV(invcf, outtsv, caller):
    '''
    Convert a VCF to TSV
    '''
    fout = open(outtsv, "w")
    vcf_in = vcf.Reader(open(invcf, "r"))
    var_ct = 0
    if caller == "gatk" or caller == "varscan":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AD", "GT"]
    elif caller == "freebayes":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AO", "GT"]
    fout.write("%s\n" % "\t".join(fields_to_extract))
    for var in vcf_in:

        chrom_ = var.CHROM
        pos_ = str(var.POS)
        ref_ = var.REF
        qual_ = str(var.QUAL)
        if caller == "varscan":
            qual_ = str(var.samples[0]["GQ"])

        dp_ = "0"
        ad_or_ao = "NA"
        ad_str = "NA"
        gt_ = "NA"
        alt_str = ""

        if var.samples[0].called == False:
            continue

        # usually need to bypass indels, however,
        # homozygous REF is considered indel by pyvcf... WTF?
        if var.is_monomorphic:
            alt_str = "."
            gt_ = "%s/%s" % (ref_, ref_)
            if caller == "freebayes":
                dp_ = var.samples[0].data.DP
                ro_ = var.samples[0].data.RO
                ad_str = str(dp_ - ro_)
            elif caller == "gatk":
                dp_ = var.samples[0].data.DP
                ad_str = "0"
            elif caller == "varscan":
                dp_ = var.INFO["ADP"]
            dp_ = str(dp_)
        else:
            alt_ = var.ALT[0]
            alt_str = "."
            if alt_ != None:
                alt_str = alt_.sequence
            if caller == "freebayes" or caller == "gatk":
                dp_ = str(var.INFO["DP"])
            else:
                dp_ = str(var.INFO["ADP"])

            gt_ = var.samples[0].gt_bases

            if caller == "gatk":
                ad_ = var.samples[0]["AD"]
                for a_ in ad_:
                    ad_str += ",%d" % a_
                ad_str = ad_str[1:]
            elif caller == "varscan":
                ad_str = str(var.samples[0]["AD"])
            else:
                ad_str = str(var.samples[0]["AO"])

        if ad_str == "NA":
            ad_str = "0"

        data_bits = [chrom_, pos_, ref_, alt_str, qual_, dp_, ad_str, gt_]
        fout.write("%s\n" % "\t".join(data_bits))
        var_ct += 1
    fout.close()
    return var_ct

def writeSampleComparisonReport(sample1, sample2, config, results):
    # unpack values used in this function
    total_compared = results['total_compared']
    ct_common = results['ct_common']
    frac_common = results['frac_common']
    comm_hom_ct = results['comm_hom_ct']
    comm_het_ct = results['comm_het_ct']
    ct_diff = results['ct_diff']
    diff_hom_ct = results['diff_hom_ct']
    diff_het_ct = results['diff_het_ct']
    diff_het_hom_ct = results['diff_het_hom_ct']
    diff_hom_het_ct = results['diff_hom_het_ct']
    diff_1sub2_ct = results['diff_1sub2_ct']
    diff_2sub1_ct = results['diff_2sub1_ct']
    judgement = results['judgement']
    short_judgement = results['short_judgement']
    # so pad numeric string to 6 spaces
    diff_hom = ("%d" % diff_hom_ct).rjust(5)
    diff_het = ("%d" % diff_het_ct).rjust(5)
    diff_het_hom = ("%d" % diff_het_hom_ct).rjust(5)
    diff_hom_het = ("%d" % diff_hom_het_ct).rjust(5)
    diff_1sub2 = ("%d" % diff_1sub2_ct).rjust(5)
    diff_2sub1 = ("%d" % diff_2sub1_ct).rjust(5)
    std_report_str = """sample1:\t%s
sample2:\t%s
variants:\t%s
depth threshold: %d
________________________________________

Positions with same genotype:   %d
breakdown:    hom: %d
          het: %d
________________________________________

Positions with diff genotype:   %d
     breakdown:
                         SAMPLE 1
                  | het  | hom  | subset
            -------+------+------+------
            het   |%s |%s |%s |
            -------+------+------+------
SAMPLE 2    hom   |%s |%s |   -  |
            -------+------+------+------
            subset|%s |   -  |   -  |
________________________________________

Total sites compared: %d
Fraction of common: %f (%d/%d)
________________________________________
CONCLUSION:
%s"""  % (sample1, sample2, config.vcf, config.dp_threshold,
          ct_common, comm_hom_ct, comm_het_ct, ct_diff, diff_het, diff_hom_het,
          diff_1sub2, diff_het_hom, diff_hom, diff_2sub1, total_compared,
          frac_common, ct_common, total_compared, judgement)
    short_report_str = """# sample1\t sample2\t DP_thresh\t FracCommon\t Same\t Same_hom\t
            Same_het\t Different\t 1het-2het\t 1het-2hom\t 1het-2sub\t 1hom-2het\t 1hom-2hom\t
            1sub-2het\t Conclusion
            %s\t%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s""" % (sample1,
                   sample2, config.dp_threshold, frac_common, ct_common, comm_hom_ct, comm_het_ct,
                   ct_diff, diff_het_ct, diff_het_hom_ct, diff_2sub1_ct, diff_hom_het_ct,
                   diff_hom_ct, diff_1sub2_ct, short_judgement)

    REPORT_PATH = "%s/%s-%s" % (config.output_dir, sample1, sample2)
    with open(REPORT_PATH, "w") as fout:
        if config.short_output:
            fout.write(short_report_str)
            if config.verbose:
                print short_report_str
        else:
            fout.write(std_report_str)
            logger.info(std_report_str)

def writeSimilarityMatrix(config, sampleSet1, sampleSet2, frac_common_matrix, total_compared_matrix):
    ''' print out grand matrix '''
    REPORT_PATH = config.report_file
    logger.info("Writing report to {}".format(REPORT_PATH))
    with open(REPORT_PATH, "w") as fout:
        for sample1 in sampleSet1:
            fout.write("\t" + sample1["name"])
        fout.write("\n")
        for sample2 in sampleSet2:
            fout.write(sample2["name"])
            for sample1 in sampleSet1:
                s = '%.4f' % frac_common_matrix[sample1["name"]][sample2["name"]]
                fout.write("\t" + s)
            fout.write("\n")
    
if __name__ == "__main__":
    main(sys.argv[1:])
