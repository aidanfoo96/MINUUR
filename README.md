<!-- ![MINUUR_Logo](assets/img/logo.png) -->

# MINUUR - Microbial INsights Using Unmapped Reads

![Main Code Base](https://img.shields.io/github/languages/top/aidanfoo96/MINUUR)
![License](https://img.shields.io/badge/license-GPLv3-blue)
![Last Commit](https://img.shields.io/github/last-commit/aidanfoo96/MINUUR)
![Open Issues](https://img.shields.io/github/issues-raw/aidanfoo96/MINUUR)
![Repo Size](https://img.shields.io/github/repo-size/aidanfoo96/MINUUR)
![MINUUR](https://github.com/aidanfoo96/MINUUR/actions/workflows/run_pipeline.yml/badge.svg)

Doi for manuscript: https://doi.org/10.12688/wellcomeopenres.19155.1

Please follow the tutorial in my Jupyter Book Available Here: https://aidanfoo96.github.io/MINUUR/ for reproduction of my analysis or to apply in your host of interest :) 

MINUUR is a snakemake pipeline I developed to extract non-host sequencing reads from mosquito whole genome sequencing data and utilise a range of metagenomic analyses to characterise potential host-associated microbes. Its application can be applied to other host-associated WGS data. MINUUR aims to leverage pre-existing WGS data to recover microbial information pertaining to host associated microbiomes.

MINUUR utilises several pieces of software in its pipeline: 
- Kraken2 to classify microbial taxa to species level from read sequences
- KrakenTools to extract classified reads pertaining to microbes for downstream analysis
- Bracken to reestimate taxonomic abundance from Kraken2 taxonomic report
- MetaPhlan3 to classify microbial taxa using marker genes
- Megahit to perform metagenome assembly 
- Quast to generate assembly statistics
- MetaBat2 to bin assembled contigs
- CheckM to assess bin quality

## Installation of Snakemake
MInUUR is run using the bioinformatics workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake is best installed using the package manager [Mamba](https://github.com/mamba-org/mamba)

Once Mamba is installed run 

`mamba create -c bioconda -c conda-forge --name snakemake snakemake`

## Installation of MINUUR
Use `git clone https://github.com/aidanfoo96/MINUUR/` and `cd MINUUR/workflow`. This is the reference point from which the pipeline will be run. See the JupyterBooks page for a full tutorial on establishing the configuration file to run this pipeline

Any feedback or bugs please open a issue or contact: aidan.foo@lstmed.ac.uk

