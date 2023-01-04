---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Introduction 

Hello! MINUUR is a metagenomics workflow I developed to pull metagenomic information from publicly available mosquito whole genome sequencing data. Data from the *Anopheles gambiae* 1000 genomes resource and various other sequencing projects on the European Nucleotide Archive produce a huge amount of short read shotgun sequencing data pertaining to the mosquito. However, many of the reads not associated to the mosquito (unmapped reads), are valuable sources of metagenomic information to us microbiome researchers. :) 

This workflow I developed serves two purposes 1. A way to reproduce a large amount of work from my PhD, which involved creating a repository of high-quality mosquito associated metagenome-assembeld genomes (MAGs) for further analysis, recovered from mosquito WGS data and 2. For those interested in following a similar approach to mine, a (hopefully) straight-forward way to perform this analysis in an organism of your choosing.

(section-label)=

## The Workflow
```{tableofcontents}
```
##### Any Questions?
Please send me an email on aidan.foo@lstmed.ac.uk.
