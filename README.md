[![Build Status](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher/statusbadges-build/icon)](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher)
[![Code Grade](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher/statusbadges-grade/icon)](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher)
[![Coverage](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher/statusbadges-coverage/icon)](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher)

# cohort-matcher #

A simple tool for determining whether two cohorts of [BAM files](https://samtools.github.io/hts-specs/SAMv1.pdf) contain reads sequenced from the same samples or patients by counting genotype matches at common SNPs. cohort-matcher is build on BAM-matcher.

BAM-matcher is most useful at comparing whole-genome-sequencing (WGS), whole-exome-sequencing (WES) and RNA-sequencing (RNA-seq) human data, but can also be customised to compare panel data or non-human data.

To compare two cohorts, run:
```
/ngs/apps/Python-2.7.8/bin/python /ngs/apps/bam-matcher/cohort_matcher.py \
        --set1 cohort1.txt --set2 cohort2.txt \
        --cache-dir `pwd`/cache --scratch-dir /scratch \
        --caller freebayes \
        --vcf /ngs/apps/bam-matcher/hg19.exome.highAF.7550.vcf \
        --reference /ngs/reference/hg19/hg19.fa \
        --freebayes-path /ngs/apps/freebayes/bin/freebayes \
        --aws /usr/bin/aws \
        --Rscript /ngs/apps/R-3.2.2/bin/Rscript \
        --samtools /ngs/apps/samtools-0.1.19/samtools \
        --output-dir output
```

which will output a series of files indicating sample similarity include:
cohort-matcher-results.txt
cohort-matcher-results.pdf
topmatches.txt
meltedResults.txt
and sample-sample comparison reports like in the cache directory:
```

BAM1:	sample1.bam 
BAM2:	sample2.bam 
depth threshold: 15 
____________________________________

Positions with same genotype:   243  
     breakdown:    hom: 51
                   het: 192
____________________________________

Positions with diff genotype:   158 
     breakdown:
                       BAM 1 
               | het  | hom  | subset 
        -------+------+------+------- 
         het   |    0 |    0 |   76  
        -------+------+------+------- 
BAM 2    hom   |    0 |    0 |   -   
        -------+------+------+------- 
         subset|   82 |   -  |   -   
____________________________________

Total sites compared: 401
Fraction of common: 0.605985 (243/401)
CONCLUSION: DIFFERENT SOURCES

```

----------------------------


# A Brief Installation Guide #

See the [wiki page](https://bitbucket.org/sacgf/bam-matcher/wiki/Installation) for detailed installation guide.



## Dependencies ##

**Python** 

(version 2.7)

**Python libraries**

* PyVCF
* ConfigParser
* Cheetah
* pysam (requires python-dev and zlib1g-dev libraries)
* fisher (requires numpy and python-dev libraries).

## Variant Callers ##

(Require at least one)

* GATK (requires Java)
* VarScan2 (requires Java and Samtools)
* Freebayes

## Installation ##

```
cd /directory/path/where/bam-matcher/is/to/be/installed/
git clone https://bitbucket.org/sacgf/bam-matcher.git
```

This provides:

```
# Python scripts and libraries
bam-matcher.py         
bammatcher_methods.py
bammatcher_exp.py
generate_example_data.py

# template files for configuration and HTML output
bam-matcher.conf.template
bam_matcher_html_template

# VCF files and example chromosome map file
1kg.exome.highAF.1511.vcf  
1kg.exome.highAF.3680.vcf
1kg.exome.highAF.7550.vcf
hg19.chromosome_map

# directory containing example BAM files
test_data/

# miscellaneous
contributors.txt
LICENSE
requirements
README.md
```

To make bam-matcher.py executable from anywhere in the system, add the directory containing bam-matcher.py to your PATH variable. e.g. add this line to your ~/.bashrc:

```
export PATH=$PATH:/path/to/bam-matcher/
```

The repository includes 3 VCF files which can be used for comparing human data (hg19/GRCh37). 

These VCF files also contain variants extracted from 1000 Genomes project which are all exonic and have high likelihood of switching between REF and ALT alleles (global allele frequency between 0.45 and 0.55). The only difference between them is the number of variants contained within.

The repository also includes several BAM files which can be used for testing (under test_data directory), as well as the expected results for various settings.


## Configuration ##

BAM-matcher requires a configuration file. The default configuration path recognised by BAM-matcher is "bam-matcher.conf" in the same directory as bam-matcher.py.

```
cd /path/to/bam-matcher/
cp bam-matcher.conf.template bam-matcher.conf
```

Then edit the file bam-matcher.conf appropriately.

If the template configuration file is missing, it can be generated by the ```--generate-config (-G)``` function. 

```
BAM-matcher.py --generate-config path_to_file_to_be_generated
```

At the very minimum, you will need to specify in the configuration file:

- ```caller:``` this is the default variant/genotype caller to use (gatk, varscan, or freebayes)
- settings for whichever caller you have chosen. For GATK, you will need to provide the path to the GATK jar file (```GATK:```); for VarScan, you will need to provide both the path to the VarScan jar file and the command to call SAMtools; for Freebayes, you will need to provide the command/path to call freebayes. GATK and VarScan will also require Java. 

- ```VCF_file:``` Specify the **FULL PATH** to the VCF file containing the variant loci to compare. Three VCF files are provided with BAM-matcher for human hg19 data.
- ```REFERENCE:``` The reference file used for mapping the reads in the input BAM files. This should also be the same version of genome reference for the VCF file.
- ```CACHE_DIR:``` You **must** supply the path to a directory with read and write permission for all users. This is used to store cached genotype data.


Most configuration settings can also be overridden at run time.

For detailed instruction on how to set up a configuration file, see the [wiki page](https://bitbucket.org/sacgf/bam-matcher/wiki/Configuration).



-----

# Running BAM-matcher #

If configured correctly, to compare two bam files, you just need to run:

```
bam-matcher.py --bam1 sample1.bam --bam2 sample2.bam -o output_report.txt
```

For detailed information on runtime arguments and parameters, see the [wiki page](https://bitbucket.org/sacgf/bam-matcher/wiki/Arguments).

See the [tutorial](https://bitbucket.org/sacgf/bam-matcher/wiki/Usage) on how to test BAM-matcher using the example data.

-----


# LICENCE #

The code is released under the Creative Commons by Attribution licence (http://creativecommons.org/licenses/by/4.0/). You are free to use and modify it for any purpose (including commercial), so long as you include appropriate attribution. 

# Citation #

*BAM-matcher: a tool for rapid NGS sample matching*

Paul P.S. Wang; Wendy T. Parker; Susan Branford; Andreas W. Schreiber
Bioinformatics 2016

[doi: 10.1093/bioinformatics/btw239](http://bioinformatics.oxfordjournals.org/content/early/2016/05/01/bioinformatics.btw239.abstract)


# Who do I talk to? #

Paul (paul.wang @ sa.gov.au)
