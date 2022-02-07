FROM ubuntu:18.04 as base

# File Author / Maintainer
LABEL org.opencontainers.image.authors="lcerdeira@gmail.com"

#configure time zone and install tzdata for base -r installation

RUN export DEBIAN_FRONTEND=noninteractive && ln -fs /usr/share/zoneinfo/Europe/London /etc/localtime

#get bits and pieces
RUN apt-get update && apt-get install --yes --no-install-recommends \
  wget \
  pkg-config \ 
  libfreetype6-dev \
  libpng-dev \
  python-matplotlib \
  locales \
  bowtie2 \
  python3-cutadapt \
  git \
  bedtools \
  cmake \
  build-essential \
  gcc-multilib \
  python3 \
  openjdk-8-jre \
  python3-pip \
  libpython3-all-dev \
  libpython2.7-dev \
  autoconf \
  automake \
  snakemake \
  make \
  clang-9 \
  g++ \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-gnutls-dev \
  libssl-dev \
  libncurses5-dev \
  ant \
  software-properties-common \
  gnupg2 \
  datamash \
  bwa \
  bzip2 \
  libgomp1  \
  gzip \
  git-lfs \
  curl \
  unzip \
  python3-setuptools \
  tzdata \
  r-base \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

WORKDIR /

#get and make megahit
RUN git clone https://github.com/voutcn/megahit.git
#RUN git submodule update --init
WORKDIR /megahit
RUN mkdir build && cd build 
RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/megahit .. && make -j4 install
RUN megahit --test && megahit --test --kmin-1pass
ENTRYPOINT ["megahit"]

#get and install QUAST
RUN wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
RUN tar -xzf quast-5.0.2.tar.gz && rm quast-5.0.2.tar.gz
WORKDIR /quast-5.0.2
RUN ./setup.py install_full

#get and install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
RUN tar -xvf samtools-1.14.tar.bz2 && rm samtools-1.14.tar.bz2
WORKDIR /samtools-1.14
RUN ./configure && make && make install

#get and install metabat2
RUN git clone https://bitbucket.org/berkeleylab/metabat.git
WORKDIR /metabat
RUN mkdir build && cd build 
RUN cmake -DCMAKE_INSTALL_PREFIX=/metabat ..  # add -DCMAKE_INSTALL_PREFIX=MY_PREFIX if needed
RUN make
RUN make install
RUN cd .. && rm -rf build

#get and install CheckM
RUN pip3 install numpy
RUN pip3 install matplotlib
RUN pip3 install pysam
RUN pip3 install checkm-genome

#get and install BEDTools (Installed via apt-get)

#get and install bwa (Installed via apt-get)

#get and install cutadapt: v1.15 (Installed via apt-get)

#get and install bowtie2: v9.3.0 (Installed via apt-get)

#get and install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip && \
  chmod 755 FastQC/fastqc && \
  ln -s $PWD/FastQC/fastqc /usr/local/bin/

#get and install Kraken2: v2.1.2
RUN wget -q http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/kraken2/
RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
RUN tar -xvf v2.1.2.tar.gz && rm v2.1.2.tar.gz
WORKDIR /Kraken2
RUN ./install_kraken2.sh Kraken2
RUN cp /Kraken2{,-build,-inspect} $HOME/bin

#get and install KrakenTools: v1.2
RUN wget https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz
RUN tar -xvf v1.2.tar.gz && rm v1.2.tar.gz
WORKDIR /KrakenTools
RUN ln -s /KrakenTools /usr/local/bin/

#get and install Metaphlan3: v3.0.13
RUN wget https://github.com/biobakery/MetaPhlAn/archive/refs/tags/3.0.14.tar.gz
RUN tar -xvf 3.0.14.tar.gz && rm 3.0.14.tar.gz
WORKDIR /Metaphlan
RUN pip install metaphlan
RUN $ metaphlan --install --bowtie2db /Metaphlan/metaphlan_databases

#get and install Bracken: v2.5.0 (note: used the v2.6.2)

RUN wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v2.6.2.tar.gz
RUN tar -xvf v2.6.2.tar.gz && rm v2.6.2.tar.gz
WORKDIR /Bracken
RUN ./install_bracken.sh Bracken
RUN cd /Bracken/src/ && make

#get and install HUMmaNn3: v3.0.0
RUN wget https://files.pythonhosted.org/packages/27/f9/d07bd76dd7dd5732c4d29d58849e96e4828c8a7dc95cf7ae58622f37591a/humann-3.0.1.tar.gz
RUN tar -xvf humann-3.0.1.tar.gz && rm humann-3.0.1.tar.gz
WORKDIR /humann
RUN pip install humann
RUN humann_databases --download chocophlan full /humann
RUN humann_databases --download uniref uniref90_diamond /humann
RUN humann_databases --download uniref uniref90_ec_filtered_diamond /humann
RUN humann_databases --download uniref uniref50_diamond /humann
RUN humann_databases --download uniref uniref50_ec_filtered_diamond /humann

#get and install R packages
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('MetBrewer')"

#get and install MInUUR
RUN git clone https://github.com/aidanfoo96/MINUUR.git
WORKDIR MINUUR

#download the DBs for minuur

WORKDIR /MINUUR/workflow/scripts
RUN ./install_db.sh

RUN export LC_ALL=C.UTF-8 && export LANG=C.UTF-8