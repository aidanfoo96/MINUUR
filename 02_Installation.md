# Installation 
In order to use this workflow, there are several things you will need to install: dependencies for the pipeline to run and several taxonomic databases! 

### Conda 
Conda is a package management system which will allow you to install and manage software environemnts. MINUUR uses conda to create isolated software environemnts to run and is also required to install the next dependency - Snakemake 

### Snakemake 
Snakemake is a workflow manager to create reproducible analyses and is the language used to write MINUUR. To install snakemake, please follow the documentation: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

# Databases 
MINUUR uses several databases to run key parts of its workflow. You are able to pick and choose what parts are relevent to your questions but the databases MINUUR requires to run a full analysis are: 
1. An indexed Kraken2 Database (alternatively a manually created Kraken2 database if you are interested in finding a specific taxa)
2. An indexed Bracken Database 
3. MetaPhlAn3 Database - you can use this in parallel or alternative to Kraken2 since MetaPhlAn3 uses marker genes instead of k-mers 
4. HUMAnN3 Database - Here a ChocoPhlAn database and UniRef90 database can be used to search for functional profiles from predicted taxa - this would be used in conjunction with MetaPhlAn3's classifications 

That's quite a few databases...

You can pick and choose what analysis to do, so feel free to omit any database you are not interested in. Certain steps can be ignored using the configuration file. See how by [configuring the workflow](04_Configuring_The_Workflow.md) 

## Installing a Database 
If the sever you are using already has one of these databases installed and ready to use, do nothing. You will specify the paths to these database in the later steps when configuring your workflow. 

If you require these databases, please go to the relevent pages and install them from there - links to each are: 

[Kraken2](https://benlangmead.github.io/aws-indexes/k2)
[Bracken](https://benlangmead.github.io/aws-indexes/k2)
[MetaPhlAn3](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0)
[Humann3](https://github.com/biobakery/humann)

