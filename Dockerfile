FROM python:2.7

# Metadata
LABEL container.base.image="python:2.7"
LABEL software.name="cohort-matcher"
LABEL software.version="0.1"
LABEL software.description="NGS Sample matching pipeline"
LABEL software.website="https://github.com/golharam/cohort-matcher"
LABEL software.documentation="README.md"
LABEL software.license="CCL 4.0"
LABEL tags="Genomics"

# System and library dependencies
RUN apt-get -y update && \
    apt-get -y install unzip build-essential cmake zlib1g-dev libncurses-dev git && \
    apt-get clean

# Get hg19
#RUN mkdir -p /ngs/reference/hg19/chromosomes
#RUN wget -O /ngs/reference/hg19/hg19.tar.gz \
#      http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz
#RUN cd /ngs/reference/hg19/chromosomes && tar zxvf ../hg19.tar.gz && \
#    cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa \
#        chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa \
#        chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chrM.fa > ../hg19.fa
#RUN rm -rf /ngs/reference/hg19/hg19.tar.gz /ngs/reference/hg19/chromosomes
# Get GRCh37
#RUN mkdir -p /ngs/reference/GRCh37
#RUN wget -O /ngs/reference/GRCh37/GRCh37.fa.gz \
#      ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && \
#    gunzip /ngs/reference/GRCh37/GRCh37.fa.gz

# Install Freebayes
RUN git config --global url.https://github.com/.insteadOf git://github.com/ && \
    git clone --recursive git://github.com/ekg/freebayes.git
RUN cd freebayes && make && make install

# Install R
## Use Debian unstable via pinning 
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list
ENV R_BASE_VERSION 3.4.1
RUN apt-get -y update && apt-get -y install r-base=${R_BASE_VERSION}*

# Samtools
RUN wget -O /samtools-1.4.tar.bz2 \
  https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 && \
  tar xvjf /samtools-1.4.tar.bz2 && rm /samtools-1.4.tar.bz2
RUN cd /samtools-1.4 && make && make install

# Other software dependencies
RUN pip install numpy awscli

# Application installation
COPY requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt && \
    pip install boto3 && \
    mkdir /cohort-matcher
COPY cohort_matcher.py /cohort-matcher/
COPY common_utils/* /cohort-matcher/common_utils/
COPY reportTopMatches.r /cohort-matcher/
COPY run_cohort_matcher.py /cohort-matcher/

# Application entry point
ENTRYPOINT ["python", "/cohort-matcher/run_cohort_matcher.py"]