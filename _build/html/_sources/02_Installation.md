# Installation 
To use MINUUR, there are several things you will need to install: dependencies for the pipeline to run (conda and Snakemake) and several taxonomic databases for microbial classification. 

### Conda 
Conda is a package management system to install and manage contained software environemnts. MINUUR uses conda to create isolated software environemnts with their specific packages and versions to run, and is also required to install the next dependency - Snakemake. If conda isn't installed on your system, please refer to the [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) to  get started.

### Snakemake 
Snakemake is a workflow manager to create reproducible analyses and was used to write MINUUR. Snakemake makes it easy to create complex and configurable workflows that should run seamlessly on other systems using isolated conda environments. Please follow the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on how to install Snakemake on your system. 

### Databases 
MINUUR uses large taxonomic databases to run several parts of its workflow. You can pick and choose which databases are relevent to your question at hand but the databases MINUUR requires to run a full analysis, or reproduce my analysis in this [paper](https://www.biorxiv.org/content/biorxiv/early/2022/08/11/2022.08.09.503283.full.pdf) are: 

1. An indexed Kraken2 Database (alternatively a manually created Kraken2 database if you are interested in finding specific taxa). Used for classifying those unmapped reads. 
2. An indexed Bracken Database. An extra step to reestimate relative abundance from your Kraken2 Classifications.
3. A MetaPhlAn3 Database - you can use this in parallel or alternatively to Kraken2 since MetaPhlAn3 uses marker genes instead of k-mers - an alternative way to classify your unmapped reads. I found Kraken2 to work better when classifying unmapped reads from mosquitoes. 

That's quite a few databases...

You can pick and choose what analysis to do and is relevent to you, so feel free to omit any database you are not interested in. Certain steps can be ignored using the configuration file (more on this later, but to skip ahead look [here](04_Configuring_The_Workflow.md)).

While obtaining all of this from external sources is a bit of a faff, the brilliant people hosting these databases will update and maintain these databases.

### Installing a Database 
If the sever you are using already has one of these databases installed and ready to use, do nothing. You will specify the paths to these database in the later steps when configuring your workflow. 

If you require these databases, please go to the relevent pages and install them from there - links to each are: 

- [Kraken2](https://benlangmead.github.io/aws-indexes/k2)
- [Bracken](https://benlangmead.github.io/aws-indexes/k2) (the Kraken and Bracken index can be used for the same analysis if downloaded from this repository)
- [MetaPhlAn3](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0)

### A Host 
MINUUR tries to recover metagenomic information from a host in question. In my case, this was from a load of mozzies. In your case...it might be something cooler. To account for different hosts, MINUUR is flexible when inputting a target host of interest. There's a caveat, however. Your host of interest must be indexed using bowtie2. Follow the instructions [here](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer). If you make your way to the top of the link provided, you will find some ready-made indexes of some common organisms (maybe yours is there!). If not, you'll need to index a chromosome-assembled fasta file of your organism in question and specify a path to your indexed reference in the configuration file.  

### Getting the Pipeline 
When you're ready, you can acquire the pipeline using `git clone https://github.com/aidanfoo96/MINUUR/` into your chosen working directory. Once you've cloned the pipeline, navigate to `cd MINUUR/workflow`. This is the reference point where you will run the pipeline. 

Try typing `snakemake -np` into your terminal once you've navigated to `MINUUR/workflow` ... you'll praobbly get an error. That's because we need to specify some data to analyse! On to the next section!

