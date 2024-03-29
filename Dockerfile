FROM ubuntu:bionic

# File Author / Maintainer
LABEL base.image="ubuntu:bionic"
LABEL dockerfile.version="1"
LABEL software="MINUUR"
LABEL software.version="0.1"
LABEL description="Pipeline to pull microbial reads from WGS data and perform metagenomic analysis"
LABEL website="https://github.com/aidanfoo96/MINUUR"
LABEL license="https://github.com/aidanfoo96/MINUUR/blob/master/LICENSE"
LABEL maintainer="Louise Cerdeira"
LABEL maintainer.email="lcerdeira@gmail.com"

# So apt doesn't ask questions during install
ARG DEBIAN_FRONTEND=noninteractive
ARG QUAST_VER="5.0.2"
ARG SAMTOOLSVER="1.14"
ARG BOWTIE2VER="2.4.5"
ARG BWAVER="0.7.17"
ARG FASTQC_VER="0.11.9"
ARG K2VER="2.1.2"
ARG KTVER="1.2"
ARG METAPH="3.0.3"
ARG BEDTOOLS_VER="2.30.0"
ARG BRACVER="2.6.2"
ARG HUMVER="3.0.1"

#get bits and pieces

RUN apt-get update && apt-get install -y --no-install-recommends \
  ant \
  autoconf \
  automake \
  build-essential \
  bzip2 \
  ca-certificates \
  cmake \
  cpanminus \
  curl \
  datamash \
  default-jre \
  g++ \
  gcc-multilib \
  git \
  gnupg2 \
  gzip \
  libboost-all-dev \
  libbz2-dev \
  libcurl4-gnutls-dev \
  libfreetype6-dev \
  libgomp1  \
  liblzma-dev \
  libncurses5-dev \
  libpng-dev \
  libpython2.7-dev \
  libpython3-all-dev \
  libssl-dev \
  locales \
  make \
  nano \
  openjdk-8-jre \
  perl \
  pkg-config \ 
  python \
  python-joblib \
  python-matplotlib \
  python-pip \
  python-setuptools \
  python-simplejson \
  python3 \
  python3-cutadapt \
  python3-distutils-extra \
  python3-pip \
  python3-setuptools \
  rsync \
  snakemake \
  software-properties-common \
  tzdata \
  unzip \
  vim \
  wget \
  zlib1g-dev \
  r-base && \
  apt-get autoclean && \
  rm -rf /var/lib/apt/lists/*

ADD https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh /tmp/miniconda.sh

RUN set -e \
  && ln -sf bash /bin/sh

RUN set -e \
  && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
  && /opt/conda/bin/conda update -n base -c defaults conda \
  && /opt/conda/bin/conda clean -ya \
  && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
  && echo '. /opt/conda/etc/profile.d/conda.sh' >> /etc/profile \
  && echo 'conda activate base' >> /etc/profile \
  && rm -f /tmp/miniconda.sh

ENV PATH /opt/conda/bin:${PATH}

RUN set -e \
  && conda config --add channels defaults \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge \
  && conda clean -ya \
  && rm -rf /root/.cache/pip

ENTRYPOINT ["/opt/conda/bin/conda"]

# RUN apt-get update && apt-get install -y software-properties-common && \
#   add-apt-repository ppa:deadsnakes/ppa && \
#   apt-get update && apt-get install -y --no-install-recommends --no-install-suggests \
#   python3.7 && \
#   python3.7 -m pip install --upgrade pip && \
#   python3.7 -m pip install numpy Cython six --force-reinstall

# FROM ubuntu:bionic AS builder_metaphlan

# # multstage build
# # labels associated with final docker image are further down
# # label the intermediate image so we can delete later 
# LABEL stage=builder_metaphlan

# # install python (>3.6) and other dependencies
# # R necessary if user wants to run unifrac function in metaphlan

# RUN apt-get update && apt-get install -y software-properties-common && \
#   add-apt-repository ppa:deadsnakes/ppa && \
#   apt-get update && apt-get install -y --no-install-recommends --no-install-suggests \
#   gcc \
#   wget \
#   python3.7 \
#   python3.7-dev \
#   python3-distutils \
#   python3-setuptools \
#   python3-pip \
#   unzip \
#   r-base=3.4.4-1ubuntu1 && \ 
#   python3.7 -m pip install pip --force-reinstall && \
#   python3.7 -m pip install numpy Cython six --force-reinstall && \
#   ln -s /usr/bin/python3.7 /usr/bin/python 

# # bowtie 2 dependency 
# RUN mkdir /usr/bin/bowtie2 && \
#   cd /usr/bin/bowtie2 && \  
#   wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download && \
#   unzip download && \
#   rm download 

# ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-2.3.3.1-linux-x86_64" \
#   LC_ALL=C

# # install metaphlan 3
# RUN python3.7 -m pip install metaphlan==3.0.3

# # download metaphlan database
# RUN metaphlan --install && \ 
#   mkdir /data 

# # build onto first stage
# # this will make the final docker image ~0.5 GB smaller 
# # after the build completes, recommend deleting the intermediate image generated from builder stage
# # can do this with the flag label=stage=filter:
# # docker image ls --filter "label=stage=builder_metaphlan"
# # docker image prune --filter "label=stage=builder_metaphlan"

# # copy over necessary bin and packages 

# COPY --from=builder_metaphlan /usr/local/bin /usr/local/bin
# #COPY --from=builder_metaphlan /usr/bin/bowtie2/ /usr/bin/bowtie2/
# COPY --from=builder_metaphlan /usr/bin/python3.7 /usr/bin/python3.7
# COPY --from=builder_metaphlan /usr/lib/python3.7 /usr/lib/python3.7
# COPY --from=builder_metaphlan /usr/lib/x86_64-linux-gnu/libexpat.so /usr/lib/x86_64-linux-gnu/libexpat.so
# COPY --from=builder_metaphlan /lib/x86_64-linux-gnu/ /lib/x86_64-linux-gnu/
# COPY --from=builder_metaphlan /usr/local/lib/python3.7/dist-packages/ /usr/local/lib/python3.7/dist-packages/
# COPY --from=builder_metaphlan /usr/lib/python3/dist-packages/* /usr/lib/python3.7/dist-packages/
# COPY --from=builder_metaphlan /usr/bin/R /usr/bin/R
# COPY --from=builder_metaphlan /usr/bin/Rscript /usr/bin/Rscript
# COPY --from=builder_metaphlan /usr/lib/R /usr/lib/R
# COPY --from=builder_metaphlan /usr/local/lib/R /usr/local/lib/R
# COPY --from=builder_metaphlan /etc/R /etc/R
# COPY --from=builder_metaphlan /usr/lib/libR.so /usr/lib/libR.so

# #ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-2.3.3.1-linux-x86_64" \ 
# #  LC_ALL=C

# # link to python, soft link doesn't get copied from intermed stage
# RUN  ln -s /usr/bin/python3.7 /usr/bin/python 

# COPY --from=build_condaforge /usr/bin/mamba /usr/bin/mamba
# COPY --from=build_condaforge /usr/bin/megahit /usr/bin/megahit
# COPY --from=build_condaforge /usr/bin/metaphlan /usr/bin/metaphlan
# COPY --from=build_condaforge /usr/bin/metabat2 /usr/bin/metabat2
# COPY --from=build_condaforge /usr/bin/checkm-genome /usr/bin/checkm-genome

#ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-2.3.3.1-linux-x86_64" \ 
# #  LC_ALL=C
#ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-2.3.3.1-linux-x86_64" \ 
# #  LC_ALL=C
#ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-2.3.3.1-linux-x86_64" \ 
# #  LC_ALL=C
#ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-2.3.3.1-linux-x86_64" \ 
# #  LC_ALL=C

# for singularity compatibility
ENV LC_ALL=C

#quast 
RUN mkdir /usr/bin/quast && \
  cd /usr/bin/quast && \
  wget https://github.com/ablab/quast/releases/download/quast_${QUAST_VER}/quast-${QUAST_VER}.tar.gz && \
  tar -xzf quast-${QUAST_VER}.tar.gz && \
  rm -rf quast-${QUAST_VER}.tar.gz && \
  cd quast-${QUAST_VER} && \
  python setup.py install

ENV PATH="$PATH:/usr/bin/quast/quast-${QUAST_VER}" \
  LC_ALL=C

#Bowtie2
RUN mkdir /usr/bin/bowtie2 && \
  cd /usr/bin/bowtie2 && \
  wget -q -O bowtie2.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2VER}/bowtie2-${BOWTIE2VER}-linux-x86_64.zip/download; \
  unzip bowtie2.zip && \
  rm bowtie2.zip

ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-${BOWTIE2VER}-linux-x86_64" \
  LC_ALL=C

#Samtools
RUN mkdir /usr/bin/samtools && \
  cd /usr/bin/samtools && \
  wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
  tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
  rm samtools-${SAMTOOLSVER}.tar.bz2 && \
  cd samtools-${SAMTOOLSVER} && \
  ./configure && \
  make

ENV PATH="$PATH:/usr/bin/samtools/samtools-${SAMTOOLSVER}" \
  LC_ALL=C

#Bwa
RUN mkdir /usr/bin/bwa &&\
  cd /usr/bin/bwa && \
  wget https://github.com/lh3/bwa/releases/download/v${BWAVER}/bwa-${BWAVER}.tar.bz2 &&\
  tar -xjf bwa-${BWAVER}.tar.bz2 &&\
  rm bwa-${BWAVER}.tar.bz2 &&\
  cd bwa-${BWAVER} &&\
  make

ENV PATH="$PATH:/usr/bin/bwa/bwa-${BWAVER}" \
  LC_ALL=C

#Bedtools
RUN cd /usr/local/bin &&\
  wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VER}/bedtools.static.binary && \
  mv bedtools.static.binary bedtools && \
  chmod +x bedtools

ENV PATH="$PATH:/usr/bin/bedtools" \
  LC_ALL=C

#fastqc
RUN mkdir /usr/bin/fastqc && \
  cd /usr/bin/fastqc && \
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip && \
  unzip fastqc_v${FASTQC_VER}.zip && \
  rm fastqc_v${FASTQC_VER}.zip && \
  chmod +x /usr/bin/fastqc/FastQC/fastqc

ENV PATH="$PATH:/usr/bin/fastqc/FastQC" \
  LC_ALL=C

#perl module required for kraken2-build
RUN cpanm Getopt::Std

#Kraken2
RUN mkdir /usr/bin/kraken2 && \
  cd /usr/bin/kraken2 && \
  wget https://github.com/DerrickWood/kraken2/archive/v${K2VER}.tar.gz && \
  tar -xzf v${K2VER}.tar.gz && \
  rm -rf v${K2VER}.tar.gz && \
  cd kraken2-${K2VER} && \
  bash install_kraken2.sh .

ENV PATH="$PATH:/usr/bin/kraken2/kraken2-${K2VER}" \
  LC_ALL=C

##### NO DATABASE INCLUDED WITH THIS DOCKER IMAGE #####
## User will need to mount a directory from their host machine that contains kraken2 database files 
## to a directory in the container (/kraken2-db exists for this purpose, but feel free to use another location)
## DL MiniKraken2_8GB database. Built from RefSeq bacteria, archaea, viral, and human libraries.
## --strip-components=1 used so that the *.k2d files end up inside /kraken2-db and not another directory
# RUN mkdir /usr/bin/kraken2/kraken2-${K2VER}/kraken2-db && \
#   cd /usr/bin/kraken2/kraken2-${K2VER}/kraken2-db && \
#   wget --no-check-certificate https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz && \
#   tar -zxf --strip-components=1 minikraken2_v2_8GB_201904.tgz && \
#   rm -rf minikraken2_v2_8GB_201904.tgz

#KrakenTools
RUN mkdir /usr/bin/krakentools && \
  cd /usr/bin/krakentools && \
  wget https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v${KTVER}.tar.gz && \
  tar -xvf v${KTVER}.tar.gz && \
  rm v${KTVER}.tar.gz

ENV PATH="$PATH:/usr/bin/krakentools/KrakenTools-${KTVER}" \
  LC_ALL=C

# install metaphlan 3
#RUN python3.7 -m pip install metaphlan==${METAPH}
#RUN python3 -m pip install metaphlan==${METAPH}

#download metaphlan database
#metaphlan --install --bowtie2db /Metaphlan/metaphlan_databases

#RUN metaphlan --install

# build onto first stage
# this will make the final docker image ~0.5 GB smaller 
# after the build completes, recommend deleting the intermediate image generated from builder stage
# can do this with the flag label=stage=filter:
# docker image ls --filter "label=stage=builder_metaphlan"
# docker image prune --filter "label=stage=builder_metaphlan"

# copy over necessary bin and packages 

# COPY --from=builder_metaphlan /usr/local/bin /usr/local/bin
# COPY --from=builder_metaphlan /usr/bin/bowtie2/ /usr/bin/bowtie2/
# COPY --from=builder_metaphlan /usr/bin/python3.7 /usr/bin/python3.7
# COPY --from=builder_metaphlan /usr/lib/python3.7 /usr/lib/python3.7
# COPY --from=builder_metaphlan /usr/lib/x86_64-linux-gnu/libexpat.so /usr/lib/x86_64-linux-gnu/libexpat.so
# COPY --from=builder_metaphlan /lib/x86_64-linux-gnu/ /lib/x86_64-linux-gnu/
# COPY --from=builder_metaphlan /usr/local/lib/python3.7/dist-packages/ /usr/local/lib/python3.7/dist-packages/
# COPY --from=builder_metaphlan /usr/lib/python3/dist-packages/* /usr/lib/python3.7/dist-packages/
# COPY --from=builder_metaphlan /usr/bin/R /usr/bin/R
# COPY --from=builder_metaphlan /usr/bin/Rscript /usr/bin/Rscript
# COPY --from=builder_metaphlan /usr/lib/R /usr/lib/R
# COPY --from=builder_metaphlan /usr/local/lib/R /usr/local/lib/R
# COPY --from=builder_metaphlan /etc/R /etc/R
# COPY --from=builder_metaphlan /usr/lib/libR.so /usr/lib/libR.so

# ENV PATH="$PATH:/usr/bin/bowtie2/bowtie2-${BOWTIE2VER}-linux-x86_64" \
#   LC_ALL=C

# link to python, soft link doesn't get copied from intermed stage
#RUN  ln -s /usr/bin/python3.7 /usr/bin/python 

#Bracken
RUN mkdir /usr/bin/bracken && \
  cd /usr/bin/bracken && \
  wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v${BRACVER}.tar.gz && \
  tar -xvf v${BRACVER}.tar.gz && \
  rm v${BRACVER}.tar.gz && \
  cd Bracken-${BRACVER} && \
  bash install_bracken.sh . && \
  cd src/ && \
  make

ENV PATH="$PATH:/usr/bin/bracken/Bracken-${BRACVER}" \
  LC_ALL=C

#Humann3
RUN mkdir /usr/bin/humann && \
  cd /usr/bin/humann && \
  wget https://files.pythonhosted.org/packages/27/f9/d07bd76dd7dd5732c4d29d58849e96e4828c8a7dc95cf7ae58622f37591a/humann-${HUMVER}.tar.gz && \
  tar -xvf humann-${HUMVER}.tar.gz && rm humann-${HUMVER}.tar.gz && \
  pip install humann 
# humann_databases --download chocophlan full /usr/bin/humann/humann-${HUMVER}/humann && \
# humann_databases --download uniref uniref90_diamond /usr/bin/humann/humann-${HUMVER}/humann && \
# humann_databases --download uniref uniref90_ec_filtered_diamond /usr/bin/humann/humann-${HUMVER}/humann && \
# humann_databases --download uniref uniref50_diamond /usr/bin/humann/humann-${HUMVER}/humann && \
# humann_databases --download uniref uniref50_ec_filtered_diamond /usr/bin/humann/humann-${HUMVER}/humann

ENV PATH="$PATH:/usr/bin/humann/humann-${HUMVER}" \
  LC_ALL=C

#install R packages
RUN Rscript -e "install.packages('tidyverse')" && \
  Rscript -e "install.packages('ggplot2')" && \
  Rscript -e "install.packages('MetBrewer')" && \
  Rscript -e "install.packages('treemapify')"

#MINUUR
RUN mkdir /usr/bin/MINUUR && \
  cd /usr/bin/MINUUR && \
  git clone https://github.com/aidanfoo96/MINUUR.git
#cd /usr/bin/MINUUR/MINUUR/workflow/scripts && \
#bash install_db.sh .
