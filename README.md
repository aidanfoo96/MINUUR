<!-- ![MINUUR_Logo](assets/img/logo.png) -->

# MINUUR - Microbial INsight Using Unmapped Reads

![Code Count](https://img.shields.io/github/languages/count/aidanfoo96/MINUUR)
![Main Code Base](https://img.shields.io/github/languages/top/aidanfoo96/MINUUR)
![Version](https://img.shields.io/badge/version-1.0-red)
![License](https://img.shields.io/badge/license-GPLv3-blue)
![Last Commit](https://img.shields.io/github/last-commit/aidanfoo96/MINUUR)
![Open Issues](https://img.shields.io/github/issues-raw/aidanfoo96/MINUUR)
![Repo Size](https://img.shields.io/github/repo-size/aidanfoo96/MINUUR)

MINUUR is still in development - any feedback is welcome! Please contact: 248064@lstmed.ac.uk or dm me on Twitter: https://twitter.com/fooheesoon

Documentation: https://aidanfoo96.github.io/MINUUR/intro.html

MINUUR is a snakemake pipeline to extract illumina short read shotgun sequencing reads and perform a range of metagenomic analyses to characterise host-associated microbes. MINUUR was developed for my PhD project and was intended to be used for the extraction and characterisation of mosquito-associated bacterial symbionts, however, its application can be applied to other host-associated WGS data. MINUUR aims to output genus & species level read classifications, functional read classifications and assembly of metagenome assembled genomes (MAGs). 

MInUUR utilises several pieces of software in its pipeline: 
- Kraken2 to classify microbial taxa to species level from read sequences
- KrakenTools to extract classified reads pertaining to microbes for downstream analysis
- Bracken to reestimate taxonomic abundance from Kraken2
- MetaPhlan3 to classify microbial taxa using marker genes (alternative to Kraken2)
- HUMMan3 to functionally profile read sequences against a ChocoPhlan and Uniref database
- Megahit to perform metagenome assembly 
- Quast to generate assembly statistics
- MetaBat2 to bin assembled contigs
- CheckM to assess bin quality

In addition, MInUUR will produce 'tidy' data in severl points of its pipeline to more easily parse and analyse downstream :) 

![workflow_fig_crop](https://user-images.githubusercontent.com/80700844/148986377-3df4a613-fd43-49e1-9093-cc631a5a0d68.png)

## Installation of Snakemake
MInUUR is run using the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake is best installed using the package manager [Mamba](https://github.com/mamba-org/mamba)

Once Mamba is installed run 

`mamba create -c bioconda -c conda-forge --name snakemake snakemake`

## Installation of MInUUR
Use `git clone https://github.com/aidanfoo96/MINUUR/` into your chosen working directory. Once closed navigate to the working directory of MINUUR using `cd MINUUR/workflow`. This is the reference point from which the pipeline will be run. 

The configuration file, `config.yaml` is where parameters of the pipeline are modified to your liking. Please see the WIKI page for a full tutorial on establishing the configuration file.

## Data Input
Inputting your data can be done in one of two ways. 
Approach 1: 
- Automatic: place your paired fastq files into the resources directory. For automatic data input, your fastq files must be named with the following format samplename_1.fastq.gz and samplename_2.fastq.gz

Approach 2: 
- Manual: Alternatively, use the sample_table.tsv file to specify a path to your paired fastq files. 

In both instances, sample names must also be listed in the samples_list.tsv file. 

#### Host Genome Requirement
MINUUR performs its analysis on unmapped reads from whole genome sequencing data. A high quality host genome is required for raw fastq inputs in order to separate reads. Download a high quality reference database of your choosing (in fasta format) and follow [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) build tutorial here to create the index.

## Database requirements
For host removal, read classification, taxonomic abundance estimation and functional read profiling, MINUUR requires several databases to be installed on the users system. All databases can be downloaded at their following repositories

#### Kraken2 Database
Download an indexed [Kraken2 and Bracken](https://benlangmead.github.io/aws-indexes/k2) database of your choosing. The standard Kraken2 database may omit for important taxa, as such MINUUR also supports classification using a larger database of Bacterial and Archaeal sequences are available from the [struo2](https://github.com/leylabmpi/Struo2) github repository, which prodvides indexed Kraken2 databases from the [GTDB](https://gtdb.ecogenomic.org/) taxonomy database [available here](http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/). Specify this database in the configuration file.   

#### MetaPhlAn3 Database
MetaPhlAn3 requires a database file containing clade specific marker genes. Installation instructions of metaphlan are found [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0).

#### HUMAnN3 Database 
[Humann3](https://github.com/biobakery/humann) requires two databases, the ChocoPhlAn database and UniRef90 translated search databases. The choices of databases and download links are available in the github page above. 

## Running Snakemake
Once the configuration file has been configured to the users choosing (see WIKI), navigate to the `workflow` directory and run `snakemake -np` to test the pipeline will run as expected. If the user is happy all rules generate the desired output, use `snakemake --cores N --use-conda` to run the pipeline, with `N` denoting the number of cores for parrelization. If no parrelization is required, use `--cores 1`. Each rule of the pipeline will be run within individual conda environments that deploy the appropriate software where required. Run this using `snakemake --cores N --use-conda`. 

## Docker image repositories & hosting
We host all of our docker images on two different repositories and periodically sync the images between the two:

  1. Docker repo for MINUUR - https://hub.docker.com/repository/docker/lcerdeira/minuur
