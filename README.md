# MINUUR
Microbial INsight Using Unmapped Reads

# About
MINUUR is a configurable snakemake pipeline to extract short shotgun sequencing reads that are unmapped to a particular organism and utilise a range of downstream metagenomic analysis steps. MINUUR differs from other metagenomic analysis pipelines by using a wide range of approaches to widely 'scavenge' biological information for host associated microbes, and to produce 'tidy' data suitable for parsing to R or Python. MINUUR utilises several programmes to extract information: 

1. Kraken2 and MetaPhlAn to classify taxa to species level
2. Humman2 to functionally characterise read sequences and extract functions of the users choosing
3. Megahit to assemble reads to metagenome assembled genomes followed by binning and quality assurance

# Installation of Snakemake
MINUUR is run using the bioinformatics workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake is best installed using the package manager [Mamba](https://github.com/mamba-org/mamba)

Once Mamba is installed run 

`mamba create -c bioconda -c conda-forge --name snakemake snakemake`

# Installation of MINUUR
Use `git clone` this repository 

### Host Genome
MINUUR separates unmapped reads from typical host whole genome sequences. A high quality host genome is required for raw fastq inputs in order to separate reads. Download a high quality reference database of your choosing (in fasta format) and follow [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) build tutorial here to create the index.

# Database requirements
For host removal, read classification and functional classification, MINUUR requires several databases to be installed on the users system. Two options are available to download the required databases 1. MINUUR provides a script to automatically install all required databases to the `resources` directory. To install the script navigate to `MINUUR/workflow/` directory and  run `scripts/install_db.sh`. This will automatically install and compile all required databases to the resources directory. 2. The user can provide their own database files and provide the directory paths in the configuration file. Instructions to install the specific databases are provided below

### Kraken2 Database
Download the indexed [Kraken2](https://benlangmead.github.io/aws-indexes/k2) database of your choosing. The standard Kraken2 database may omit for important taxa, as such MINUUR also supports classification using a larger database of Bacterial and Archaeal sequences are available from the [struo2](https://github.com/leylabmpi/Struo2) github repository, which prodvides indexed Kraken2 databases from the [GTDB](https://gtdb.ecogenomic.org/) taxonomy database [available here](http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/).  

### Bracken Database
For reestimatation of the Kraken2 output, the used Kraken database is required to be build for Bracken. Instructions on Bracken-build are available [here](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual).

### MetaPhlAn3 Database
MetaPhlAn3 requires a database file containing clade specific marker genes. Installation instructions of metaphlan are found [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0).

### HUMAnN3 Database 
[Humann3](https://github.com/biobakery/humann) requires two databases, the ChocoPhlAn database and translated search databases. The choices of databases and download links are available on the github page above. 

