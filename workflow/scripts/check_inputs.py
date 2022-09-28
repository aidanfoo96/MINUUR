#!/usr/bin/env python3

import sys
sys.stderr = open(snakemake.log[0], "w")
import os
import pandas as pd
import numpy as np


sample_list = pd.read_csv(snakemake.params['sample_list'], sep  = "\t")
automatic_input = snakemake.params['automatic_input']

# Check sample names given in the fastq_samples.tsv file match the given fastq files

if automatic_input is False: 
    sample_table = pd.read_csv(snakemake.params['sample_table'], sep = "\t")
    assert np.isin(sample_table['sampleID'], sample_list['sampleID']).all(), f"samples in your sample table are not present in your sample list. Please ensure all samples match"
else:
    for sample in sample_list['sampleID']:
        for num in [1, 2]:
            fastq_sample = f"../resources/{sample}_{num}.fastq.gz"
            assert os.path.isfile(fastq_sample), "your sample names given in 'samples_fastq.tsv' do not match your fastq files. Please rename your files accrordingly"

print("Inputs Passed :D")