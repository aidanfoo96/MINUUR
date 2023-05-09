#!/usr/bin/env python
# coding: utf-8

# # Introduction 
# 
# Hello! MINUUR is a metagenomics workflow I developed to pull metagenomic information from publicly available mosquito whole genome sequencing data. Data from the *Anopheles gambiae* 1000 genomes resource and various other sequencing projects on the European Nucleotide Archive produce a huge amount of short read shotgun sequencing data pertaining to the mosquito. However, many of the reads not associated to the mosquito (unmapped reads), are valuable sources of metagenomic information to us microbiome researchers. :) 
# 
# This workflow I developed serves two purposes 1. A way to reproduce a large amount of work from my PhD, which involved creating a repository of high-quality mosquito associated metagenome-assembeld genomes (MAGs) for further analysis, recovered from mosquito WGS data and 2. For those interested in following a similar approach to mine, a (hopefully) straight-forward way to perform this analysis in an organism of your choosing.
# 
# One of the biggest points to consider when running this type of analysis is whether the information you are recovering pertains to "true" symbionts or artifacts from sequencing (i.e index hopping, contamination from library prep etc) - see these papers [1](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1004437) [2](https://www.nature.com/articles/s41598-022-13269-z). Making this distinction is difficult, but some steps I would suggest would be to cross reference your MAGs with symbionts previously identified within your host domain and commonly identified contaminants, and placing your genomes within the wider context of the species population. Here you could use FastANI or MASH to identify their level of similarity with one another and other host-associated microbes (if this information is available to you).
# 
# (section-label)=
# 
# ## The Workflow
# ```{tableofcontents}
# ```
# ##### Any Questions?
# Please send me an email on aidan.foo@lstmed.ac.uk.
