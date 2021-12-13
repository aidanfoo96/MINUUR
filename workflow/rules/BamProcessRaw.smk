################ Process RAW BAM FILES ################
    # Samtools: (Li et al., 2009)
    # BEDTools: (Quinlan & Hall, 2010)
#####################################################


def get_input(wildcards): 
    """
        Get sample names from sample_table (specified in Snakefile)
    """
    unit = sample_table.loc[wildcards.sample]
    if pd.isna(unit["fastq1"]):
        return "data/{bam}.bam".format(
            bam = unit["sampleID"]
        )

    else: 
        return expand(
            "data/{sample}_{{read}}.fastq.gz".format(
                sample = unit["sampleID"],
            ),
            read = ["1", "2"]
        )
 

rule bam_file_statistic: 
    """
        Generate statistics of raw bam files
    """
    output:
        txt = "../results/bam_stats/{sample}_stats.txt",
    input: 
        get_input
    shell: 
        r"""
            samtools stats {input} > {output.txt}
         """


rule print_bam_statistics: 
    """
        Extract total, mapped and unmapped reads 
    """
    output: 
        txt = "../results/alignment_stats/{sample}_align_stats.txt",
    input: 
        sample = "../results/bam_stats/{sample}_stats.txt",
    shell: 
        "sed -n '8p;14p;16p' {input.sample} > {output.txt}" 



rule concatenate_alignment_statistics: 
    """
        Concatenate total, mapped and unmapped read stats from all samples specified
    """
    output: 
        txt = "../results/alignment_stats/concatenated_alignment_statistics.txt",
    input: 
        sample = expand("../results/alignment_stats/{sample}_align_stats.txt", sample=samples),
    params: 
        filename = "FILENAME",
    shell: 
        r"""
            awk '{{print $0 "\t" {params.filename}}}' {input.sample} > {output.txt}
        
         """


rule extract_unmapped_bam:    
    """
        Extract unmapped reads from the raw bam file
    """
    output: 
        bam = "../results/unmapped_bam/{sample}_unmapped.bam",
    input: 
        get_input
    conda: 
        "../envs/bam_processing_env.yaml",
    shell: 
        "samtools view -b -f 4 {input} > {output.bam} "


rule bam_to_fastq: 
    """
        Convert the bam file to a paired fastq file 

    """
    output: 
        fastq1 = "../results/unmapped_fastq/{sample}_unmapped_1.fastq",
        fastq2 = "../results/unmapped_fastq/{sample}_unmapped_2.fastq",
    input: 
        sample = "../results/unmapped_bam/{sample}_unmapped.bam",
    conda: 
        "../envs/bam_processing_env.yaml",
    shell: 
        "bamToFastq -i {input.sample} -fq {output.fastq1} -fq2 {output.fastq2}"