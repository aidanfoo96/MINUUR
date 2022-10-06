################ Process RAW FASTQ FILES ################
    # Samtools: (Li et al., 2009)
    # BEDTools: (Quinlan & Hall, 2010)
#####################################################


#-------------------------------------------------------------------#
rule BamStatsFFQ: 
    """
        Generate Statistics of aligned bam files
    """
    output: 
        txt = "../results/bam_stats_ffq/{sample}_stats.txt",
    input: 
        sample = "../results/aligned_bam/{sample}_sorted.bam",
    log: 
        "logs/samtool_stats_ffq/{sample}.log", 
    wrapper: 
        "0.79.0/bio/samtools/stats"


#-------------------------------------------------------------------#
rule PrintBamStatsFFQ: 
    """
        Extract total, mapped and unmapped reads 
    """
    output: 
        txt = "../results/alignment_stats_ffq/{sample}_align_stats.txt",
    input: 
        sample = "../results/bam_stats_ffq/{sample}_stats.txt",
    shell: 
        "sed -n '8p;14p;16p' {input.sample} > {output.txt}" 


#-------------------------------------------------------------------#
rule ConcatenateBamStat: 
    """
        Concatenate all mapping statistic files 
    """
    output: 
        txt = "../results/alignment_stats_ffq/concatenated_alignment_statistics.txt",
    input: 
        sample = expand("../results/alignment_stats_ffq/{sample}_align_stats.txt", sample=samples),
    params: 
        filename = "FILENAME",
    shell: 
        r"""
            awk '{{print $0 "\t" {params.filename}}}' {input.sample} > {output.txt}
        
         """


#-------------------------------------------------------------------#
rule ExtractUnmappedReads:  
    """
        Extract unmapped reads from bam 
    """
    output: 
        bam = temporary("../results/unmapped_bam_ffq/{sample}_unmapped.bam"),
    input: 
        sample = "../results/aligned_bam/{sample}_sorted.bam",
    conda: 
        "../envs/bam_processing_env.yaml",
    benchmark: 
        "benchmarks/02_BamProcessFromFastq/{sample}.BAM_separation.benchmark.txt",
    log: 
        "logs/samtools_extract_ffq/{sample}.log", 
    shell: 
        "samtools view -b -f 12 -F 256 {input.sample} > {output.bam} 2> {log}"


#-------------------------------------------------------------------#
rule ConvertUnmappedReadsToFastq: 
    """
        Convert bam to paried Fastq
    """
    output: 
        fastq1 = temporary("../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq"),
        fastq2 = temporary("../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq"),
    input: 
        sample = "../results/unmapped_bam_ffq/{sample}_unmapped.bam",
    benchmark: 
        "benchmarks/02_BamProcessFromFastq/{sample}.BAMtoFastq.benchmark.txt",
    conda: 
        "../envs/bam_processing_env.yaml",
    shell: 
        "bamToFastq -i {input.sample} -fq {output.fastq1} -fq2 {output.fastq2}"


#-------------------------------------------------------------------#
rule GzipFastqFiles:
    """
        Gzip fastq_files 
    """
    output: 
        fastq1gz = "../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq.gz",
        fastq2gz = "../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq.gz",
    input:
        fastq1 = "../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq",
        fastq2 = "../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq",
    conda: 
        "../envs/bam_processing_env.yaml",
    shell: 
        r"""
            gzip -nc {input.fastq1} > {output.fastq1gz}
            gzip -nc {input.fastq2} > {output.fastq2gz}

        """
