MINUUR is still in development - any feedback is welcome! Please contact: 248064@lstmed.ac.uk or dm me on Twitter: https://twitter.com/fooheesoon

# MINUUR - Microbial INsight Using Unmapped Reads
MINUUR is a snakemake pipeline to extract unmapped whole genome shotgun sequencing reads and utilise a range of metagenomic analyses to characterise host-associated microbes. Orginally, MINUUR was intended to be used for the extraction of mosquito-associated bacterial symbionts, however, its application can be used in other host-associated WGS data. MINUUR aims to leverage pre-existing WGS data to 'scavenge' for microbial information pertaining to host associated microbiomes - the key advantage being metagenomic reads as inputs to produce species level classifications, functional inference and assembly of metagenome assembled genomes (MAGs). 

MINUUR utilises several softwares in its pipeline: 
- Kraken2 to classify microbial taxa to species level from read sequences
- KrakenTools to extract classified reads pertaining to microbes for downstream analysis
- Bracken to reestimate taxonomic abundance from Kraken2 taxonomic report
- MetaPhlan3 to classify microbial taxa using marker genes
- HUMMan3 to functionally profile read sequences against the ChocoPhlan and Uniref databases
- Megahit to perform metagenome assembly 
- Quast to generate assembly statistics
- MetaBat2 to bin assembled contigs
- CheckM to assess bin quality

In addition, MINUUR will produce 'tidy' data suitable for parsing to R or Python.

![workflow_fig_crop](https://user-images.githubusercontent.com/80700844/148986377-3df4a613-fd43-49e1-9093-cc631a5a0d68.png)

## Installation of Snakemake
MINUUR is run using the bioinformatics workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake is best installed using the package manager [Mamba](https://github.com/mamba-org/mamba)

Once Mamba is installed run 

`mamba create -c bioconda -c conda-forge --name snakemake snakemake`

## Installation of MINUUR
Use `git clone https://github.com/aidanfoo96/MINUUR/` and `cd MINUUR/workflow`. This is the reference point from which the pipeline will be run. See the WIKI page for a full tutorial on establishing the configuration file to run the pipeline

#### Host Genome
MINUUR separates unmapped reads from typical host whole genome sequences. A high quality host genome is required for raw fastq inputs in order to separate reads. Download a high quality reference database of your choosing (in fasta format) and follow [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) build tutorial here to create the index.

## Database requirements
For host removal, read classification and functional classification, MINUUR requires several databases to be installed on the users system. Two options are available to download the required databases 1. MINUUR provides a script to automatically install all required databases to the `resources` directory. To install the script navigate to `MINUUR/workflow/` directory and  run `scripts/install_db.sh`. This will automatically install and compile all required databases to the resources directory. 2. The user can provide their own database files and provide the directory paths in the configuration file. Instructions to install the specific databases are provided below

#### Kraken2 Database
Download the indexed [Kraken2](https://benlangmead.github.io/aws-indexes/k2) database of your choosing. The standard Kraken2 database may omit for important taxa, as such MINUUR also supports classification using a larger database of Bacterial and Archaeal sequences are available from the [struo2](https://github.com/leylabmpi/Struo2) github repository, which prodvides indexed Kraken2 databases from the [GTDB](https://gtdb.ecogenomic.org/) taxonomy database [available here](http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/).  

#### Bracken Database
For reestimatation of the Kraken2 output, the used Kraken database is required to be build for Bracken. Instructions on Bracken-build are available [here](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual).

#### MetaPhlAn3 Database
MetaPhlAn3 requires a database file containing clade specific marker genes. Installation instructions of metaphlan are found [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0).

#### HUMAnN3 Database 
[Humann3](https://github.com/biobakery/humann) requires two databases, the ChocoPhlAn database and translated search databases. The choices of databases and download links are available on the github page above. 

