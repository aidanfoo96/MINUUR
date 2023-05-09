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

MINUUR utilises: 
- KRAKEN2: Classify taxa from unmapped read sequences
- KrakenTools: extract classified reads for downstream analysis
- BRACKEN: reestimate taxonomic abundance from KRAKEN2
- MetaPhlan3: Classify taxa using marker genes
- MEGAHIT: Metagenome assemblies using unmapped reads
- QUAST: Assembly statistics from MEGAHIT assemblies
- MetaBat2: Bin contiguous sequences from MEGAHIT
- CheckM: Assess bin quality from MetaBat2

## Installation of Snakemake
MINUUR is run using the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Snakemake is best installed using the package manager [Mamba](https://github.com/mamba-org/mamba)

Once Mamba is installed run 

`mamba create -c bioconda -c conda-forge --name snakemake snakemake`

## Installation of MINUUR
Use `git clone https://github.com/aidanfoo96/MINUUR/` and `cd MINUUR/workflow`. This is the reference point from which the pipeline will be run. See the JupyterBooks page for a full tutorial on establishing the configuration to run this pipeline. 

### Update 09/05/2023: 
- Added Github actions 
- Dummy dataset now included in `workflow/data`, tutorial for running this is included in the JupyterBooks page. Use this to ensure the pipeline works on your machine. 
- Added the option to run BUSCO to help assess eukaryotic contamination in MAGs

Any feedback or bugs please open an issue or contact: aidan.foo@lstmed.ac.uk

