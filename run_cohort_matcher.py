from __future__ import print_function
''' This is a wrapper script for Docker '''

import multiprocessing
import os
import shlex
import subprocess
from argparse import ArgumentParser

from common_utils.s3_utils import download_file, upload_file, download_folder, upload_folder
from common_utils.job_utils import generate_working_dir, delete_working_dir, uncompress


def download_reference(s3_path, working_dir):
    """
    Downloads reference folder that has been configured to run with Isaac
    :param s3_path: S3 path that the folder resides in
    :param working_dir: working directory
    :return: local path to the folder containing the reference
    """

    reference_folder = os.path.join(working_dir, 'reference')

    try:
        os.mkdir(reference_folder)
    except Exception as e:
        pass

    download_folder(s3_path, reference_folder)

    # Update sorted reference
    update_sorted_reference(reference_folder)

    return reference_folder


def download_fastq_files(fastq1_s3_path, fastq2_s3_path, working_dir):
    """
    Downlodas the fastq files
    :param fastq1_s3_path: S3 path containing FASTQ with read1
    :param fastq2_s3_path: S3 path containing FASTQ with read2
    :param working_dir: working directory
    :return: local path to the folder containing the fastq
    """
    fastq_folder = os.path.join(working_dir, 'fastq')

    try:
        os.mkdir(fastq_folder)
    except Exception as e:
        pass

    local_fastq1_path = download_file(fastq1_s3_path, fastq_folder)
    local_fastq2_path = download_file(fastq2_s3_path, fastq_folder)

    # Isaac requires the fastqs to be symlinked as lane1_read1.fastq.gz and lane1_read2.fastq.gz
    os.symlink(local_fastq1_path, os.path.join(fastq_folder, 'lane1_read1.fastq.gz'))
    os.symlink(local_fastq2_path, os.path.join(fastq_folder, 'lane1_read2.fastq.gz'))

    return fastq_folder


def upload_bam(bam_s3_path, local_folder_path):
    """
    Uploads results folder containing the bam file (and associated output)
    :param bam_s3_path: S3 path to upload the alignment results to
    :param local_folder_path: local path containing the alignment results
    """

    upload_folder(bam_s3_path, local_folder_path)

def run_cohort_matcher(bam_sheet1, bam_sheet2, reference1, reference2, working_dir, output_prefix):
    """
    Runs Cohort-matcher
    :param bam_sheet1: local path to bam sheet1
    :param bam_sheet2: local path to bam sheet2
    :param reference1: hg19 or GRCh37
    :param reference2: hg19 or GRCh37
    :param working_dir: working directory
    :param output_prefix: output prefix
    :return: path to results
    """

    os.chdir(working_dir)
    cache_dir = os.path.join(working_dir, 'cache')
    try:
        os.mkdir(cache_dir)
    except Exception as e:
        pass
    
    max_jobs = multiprocessing.cpu_count()
    if reference1 == 'hg19':
        ref = os.path.join(working_dir, 'hg19.fa')
        vcf = os.path.join(working_dir, 'hg19.exome.highAF.7550.vcf')
        if reference2 == 'hg19':
            ref2 = ''
        else:
            ref2 = '-R2 %s -V2 %s -CM %s' % (os.path.join(working_dir, 'GRCh37.fa'),
                                             os.path.join(working_dir, '1kg.exome.highAF.7550.vcf'),
                                             os.path.join(working_dir, 'hg19.chromosome_map'))
    else:
        ref = os.path.join(working_dir, 'GRCh37.fa')
        vcf = os.path.join(working_dir, '1kg.exome.highAF.7550.vcf')
        if reference2 == 'GRCh37':
            ref2 = ''
        else:
            ref2 = '-R2 %s -V2 %s -CM %s' % (os.path.join(working_dir, 'hg19.fa'),
                                             os.path.join(working_dir, 'hg19.exome.highAF.7550.vcf'),
                                             os.path.join(working_dir, 'hg19.chromosome_map'))

    cmd = '/cohort_matcher.py --set1 %s --set2 %s --cache-dir %s --scratch-dir %s ' \
          '--caller freebayes --max-jobs %d -R %s -V %s %s ' \
          '--freebayes-path /usr/local/bin/freebayes --aws /usr/local/bin/aws ' \
          '--Rscript /usr/bin/Rscript --samtools /usr/local/bin/samtools --output_prefix %s' % \
          (bam_sheet1, bam_sheet2, cache_dir, working_dir, max_jobs, ref, vcf, ref2, output_prefix)
    print ("Running: %s" % cmd)
    subprocess.check_call(shlex.split(cmd))
    return working_dir

def parseArguments():
    argparser = ArgumentParser()

    file_path_group = argparser.add_argument_group(title='File paths')
    file_path_group.add_argument('--set1_s3_path', type=str, help='S3 path to first set of samples', required=True)
    file_path_group.add_argument('--set2_s3_path', type=str, help='S3 path to second set of samples', required=True)
    file_path_group.add_argument('--set1_reference', type=str, help='Reference genome for set1', required=True,
				 choices=['hg19', 'GRCh37'])
    file_path_group.add_argument('--set2_reference', type=str, help='Reference genome for set2', required=True,
                                 choices=['hg19', 'GRCh37'])
    file_path_group.add_argument('--s3_output_folder_path', type=str, help='S3 path for output files', required=True)

    run_group = argparser.add_argument_group(title='Run command args')
    run_group.add_argument('--output_prefix', type=str, help='Output prefix')

    argparser.add_argument('--working_dir', type=str, default='/scratch', help="Working directory (default: /scratch)")

    return argparser.parse_args()

def main():
    args = parseArguments()

    working_dir = generate_working_dir(args.working_dir)

    # Download fastq files and reference files
    print ('Downloading bam sheets')
    set1_bamsheet = download_file(args.set1_s3_path, working_dir)
    set2_bamsheet = download_file(args.set2_s3_path, working_dir)

    # Download reference bundles
    if args.set1_reference == 'hg19' or args.set2_reference == 'hg19':
        download_file('s3://bmsrd-ngs-repo/reference/hg19-cohort-matcher.tar.bz2', working_dir)
        uncompress(os.path.join(working_dir, 'hg19-cohort-matcher.tar.bz2'), working_dir)
    if args.set2_reference == 'GRCh37' or args.set2_reference == 'GRCh37':
        download_file('s3://bmsrd-ngs-repo/reference/GRCh37-cohort-matcher.tar.bz2', working_dir)
        uncompress(os.path.join(working_dir, 'GRCh37-cohort-matcher.tar.bz2', working_dir))

    # Run cohort-matcher
    print ('Running cohort-matcher')
    output_folder_path = run_cohort_matcher(set1_bamsheet, set2_bamsheet, args.set1_reference,
                                            args.set2_reference, working_dir, args.output_prefix)
    print ('Uploading results to %s' % args.output_s3_folder_path)
    #upload_bam(args.bam_s3_folder_path, bam_folder_path)
    print('Cleaning up working dir')
    delete_working_dir(working_dir)
    print ('Completed')

if __name__ == '__main__':
    main()