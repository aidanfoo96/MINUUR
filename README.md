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

# Database requirements
For host removal, read classification and functional classification, MINUUR requires several databases to be installed on the users system

## Bowtie2 Indexed Host Database
Download a high quality host genome of choosing. Follow [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) build tutorial here to create the index

## Kraken2 Database
Download the indexed [Kraken2](https://benlangmead.github.io/aws-indexes/k2) database of your choosing. The standard Kraken2 database may omit for important taxa, as such MINUUR also supports classification using a larger database of Bacterial and Archaeal sequences are available from the [struo](https://github.com/leylabmpi/Struo) github repository, which prodvides indexed Kraken2 databases from the [GTDB](https://gtdb.ecogenomic.org/) taxonomy database. 


# Tutorial
