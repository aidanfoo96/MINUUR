################ Process RAW FASTQ FILES ################
    # Samtools: (Li et al., 2009)
    # BEDTools: (Quinlan & Hall, 2010)
#####################################################


rule bam_file_statsistic_fq: 
    """
        Generate Statistics of aligned bam files
    """
    output: 
        txt = "../results/bam_stats_ffq/{sample}_stats.txt",
    input: 
        sample = "../results/aligned_bam/{sample}_sorted.bam",
    wrapper: 
        "0.79.0/bio/samtools/stats"


rule print_bam_statistics_fq: 
    """
        Extract total, mapped and unmapped reads 
    """
    output: 
        txt = "../results/alignment_stats_ffq/{sample}_align_stats.txt",
    input: 
        sample = "../results/bam_stats_ffq/{sample}_stats.txt",
    shell: 
        "sed -n '8p;14p;16p' {input.sample} > {output.txt}" 


rule concatenate_alignment_statistics_fq: 
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


rule extract_unmapped_bam_from_fastq:  
    """
        Extract unmapped reads from bam 
    """
    output: 
        bam = "../results/unmapped_bam_ffq/{sample}_unmapped.bam",
    input: 
        sample = "../results/aligned_bam/{sample}_sorted.bam",
    conda: 
        "../envs/bam_processing_env.yaml",
    shell: 
        "samtools view -b -f 12 -F 256 {input.sample} > {output.bam} "



rule bam_to_fastq_ffq: 
    """
        Convert bam to paried Fastq
    """
    output: 
        fastq1 = "../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq",
        fastq2 = "../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq",
    input: 
        sample = "../results/unmapped_bam_ffq/{sample}_unmapped.bam",
    conda: 
        "../envs/bam_processing_env.yaml",
    shell: 
        "bamToFastq -i {input.sample} -fq {output.fastq1} -fq2 {output.fastq2}"


rule count_reads_fastq_ffq: 
    """
        Count number of reads from fastq file
    """
    output: 
        fq_count = "../results/alignment_stats_ffq/{sample}_count.txt",
        fq_count2 = "../results/alignment_stats_ffq/{sample}_count2.txt",

    input: 
        fastq = expand("../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq", sample = samples),
        fastq2 = expand("../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq", sample = samples),

    shell: 
        r"""
            echo $(( $(wc -l <{input.fastq}) / 4 )) > {output.fq_count} 
            echo $(( $(wc -l <{input.fastq2}) / 4 )) > {output.fq_count2}

         """