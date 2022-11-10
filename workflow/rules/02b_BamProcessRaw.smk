################ Process RAW BAM FILES ################
    # Samtools: (Li et al., 2009)
    # BEDTools: (Quinlan & Hall, 2010)
#####################################################


#-------------------------------------------------------------------#
def GetBamInput(wildcards):
    """
    Function to determine if input is raw fastq or BAM
    """
    if config["input_config"]["automatic"]:
        units = pd.read_csv(config["samples"], sep = "\t")
        units = (
            units.assign(BamPath=f"../resources/" + units["sampleID"] + ".bam")
            .set_index("sampleID")
        )
    else:
        assert os.path.isfile(
            config["input_config"]["manual"]
            ), f"config['input_config']['manual'] doesn't exist. Please create one or use the 'automatic option' with the correct naming scheme specified in the WIKI page"
        units = pd.read_csv(config["input_config"]["manual"], sep = "\t", index_col="sampleID")
    
    u = units.loc[wildcards.sample, ["BamPath"]].dropna()

    return [f"{u.BamPath}"]
 

#-------------------------------------------------------------------#
rule bam_file_statistic: 
    """
        Generate statistics of raw bam files
    """
    output:
        txt = "../results/bam_stats/{sample}_stats.txt",
    input: 
        GetBamInput
    log: 
        "logs/samtool_stats/{sample}.log", 
    shell: 
        r"""
            samtools stats {input} > {output.txt} 2> {log}
         """


#-------------------------------------------------------------------#
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


#-------------------------------------------------------------------#
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


#-------------------------------------------------------------------#
rule extract_unmapped_bam:    
    """
        Extract unmapped reads from the raw bam file
    """
    output: 
        bam = "../results/unmapped_bam/{sample}_unmapped.bam",
    input: 
        GetBamInput
    conda: 
        "../envs/bam_processing_env.yaml",
    log: 
        "logs/samtools_extract/{sample}.log", 
    shell: 
        "samtools view -b -f 4 {input} > {output.bam} 2> {log}"


#-------------------------------------------------------------------#
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