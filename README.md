<!-- ![MINUUR_Logo](assets/img/logo.png) -->

# MInUUR - Microbial INsight Using Unmapped Reads

![Code Count](https://img.shields.io/github/languages/count/aidanfoo96/MINUUR)
![Main Code Base](https://img.shields.io/github/languages/top/aidanfoo96/MINUUR)
![Version](https://img.shields.io/badge/version-1.0-red)
![License](https://img.shields.io/badge/license-GPLv3-blue)
![Last Commit](https://img.shields.io/github/last-commit/aidanfoo96/MINUUR)
![Open Issues](https://img.shields.io/github/issues-raw/aidanfoo96/MINUURT)
![Repo Size](https://img.shields.io/github/repo-size/aidanfoo96/MINUUR)

MInUUR is still in development - any feedback is welcome! Please contact: 248064@lstmed.ac.uk or dm me on Twitter: https://twitter.com/fooheesoon

MInUUR is a snakemake pipeline to extract unmapped whole genome shotgun sequencing reads and utilise a range of metagenomic analyses to characterise host-associated microbes. Orginally, MInUUR was intended to be used for the extraction of mosquito-associated bacterial symbionts, however, its application can be applied to other host-associated WGS data. MInUUR aims to leverage pre-existing WGS data to 'scavenge' for microbial information pertaining to host associated microbiomes - the key advantage being metagenomic reads as inputs to produce genus & species level classifications, functional inference and assembly of metagenome assembled genomes (MAGs). 

MInUUR utilises several pieces of software in its pipeline: 
- Kraken2 to classify microbial taxa to species level from read sequences
- KrakenTools to extract classified reads pertaining to microbes for downstream analysis
- Bracken to reestimate taxonomic abundance from Kraken2 taxonomic report
- MetaPhlan3 to classify microbial taxa using marker genes
- HUMMan3 to functionally profile read sequences against the ChocoPhlan and Uniref databases
- Megahit to perform metagenome assembly 
- Quast to generate assembly statistics
- MetaBat2 to bin assembled contigs
- CheckM to assess bin quality

In addition, MInUUR will produce 'tidy' data suitable for parsing to R or Python.

![workflow_fig_crop](https://user-images.githubusercontent.com/80700844/148986377-3df4a613-fd43-49e1-9093-cc631a5a0d68.png)

## Installation of Snakemake
MInUUR is run using the bioinformatics workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake is best installed using the package manager [Mamba](https://github.com/mamba-org/mamba)

Once Mamba is installed run 

`mamba create -c bioconda -c conda-forge --name snakemake snakemake`

## Installation of MInUUR
Use `git clone https://github.com/aidanfoo96/MINUUR/` and `cd MINUUR/workflow`. This is the reference point from which the pipeline will be run. See the WIKI page for a full tutorial on establishing the configuration file to run the pipeline

#### Host Genome
MINUUR separates unmapped reads from typical host whole genome sequences. A high quality host genome is required for raw fastq inputs in order to separate reads. Download a high quality reference database of your choosing (in fasta format) and follow [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) build tutorial here to create the index.

## Database requirements
For host removal, read classification, taxonomic abundance estimation and functional read profiling, MINUUR requires several databases to be installed on the users system. To recreate the analysis of our paper, all databases can be downloaded at their following repositories

#### Kraken2 Database
Download the indexed [Kraken2](https://benlangmead.github.io/aws-indexes/k2) database of your choosing. The standard Kraken2 database may omit for important taxa, as such MINUUR also supports classification using a larger database of Bacterial and Archaeal sequences are available from the [struo2](https://github.com/leylabmpi/Struo2) github repository, which prodvides indexed Kraken2 databases from the [GTDB](https://gtdb.ecogenomic.org/) taxonomy database [available here](http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/).  

#### Bracken Database
For reestimatation of the Kraken2 output, the used Kraken database is required to be build for Bracken. Instructions on Bracken-build are available [here](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual).

#### MetaPhlAn3 Database
MetaPhlAn3 requires a database file containing clade specific marker genes. Installation instructions of metaphlan are found [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0).

#### HUMAnN3 Database 
[Humann3](https://github.com/biobakery/humann) requires two databases, the ChocoPhlAn database and translated search databases. The choices of databases and download links are available on the github page above. 

## Running Snakemake
Once the configuration file has been configured to the users choosing (see WIKI), navigate to the `workflow` directory and run `snakemake -np` to test the pipeline will run as expected. If the user is happy all rules generate the desired output, use `snakemake --cores N --use-conda` to run the pipeline, with `N` denoting the number of cores for parrelization. If no parrelization is required, use `--cores 1`. Each rule of the snakemake pipeline can be run within individual conda environments that deploy the software when required. Run this using `snakemake --cores N --use-conda`

## Docker image repositories & hosting
We host all of our docker images on two different repositories and periodically sync the images between the two:

  1. Docker repo for MINUUR - https://hub.docker.com/repository/docker/lcerdeira/minuur
