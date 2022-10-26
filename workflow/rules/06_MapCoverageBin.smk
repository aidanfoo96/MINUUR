################ Map Reads and Bin ##############################
    # bwa: (Li & Durbin, 2009)
    # samtools: (Li et al., 2009)
    # metabat2: (Kang et al., 2019)
##################################################################


#-------------------------------------------------------------------#
def megahit_input_1(wildcards): 
    if config["MetagenomeAssm"]["UseKrakenExtracted"]:             
        return "../results/kraken_taxon_extract/{sample}_extract_1.fastq"
    elif config["input_type"]["fastq"]:
        return "../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq"
    elif config["input_type"]["bam"]: 
        return "../results/unmapped_fastq/{sample}_unmapped_1.fastq"
    

#-------------------------------------------------------------------#
def megahit_input_2(wildcards): 
    if config["MetagenomeAssm"]["UseKrakenExtracted"]: 
        return "../results/kraken_taxon_extract/{sample}_extract_2.fastq"
    elif config["input_type"]["fastq"]:
        return "../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq"
    elif config["input_type"]["bam"]: 
        return "../results/unmapped_fastq/{sample}_unmapped_2.fastq"


#-------------------------------------------------------------------#
rule IndexMegahitContigs:
    """
        Index megahit assembled contigs using bwa
    """
    output: 
        "../results/megahit_assm/{sample}_assm/{sample}.amb",
        "../results/megahit_assm/{sample}_assm/{sample}.ann",
        "../results/megahit_assm/{sample}_assm/{sample}.bwt",
        "../results/megahit_assm/{sample}_assm/{sample}.pac",
        "../results/megahit_assm/{sample}_assm/{sample}.sa",
    input: 
        "../results/megahit_assm/{sample}_assm/final.contigs.fa",
    log:
        "logs/bwa_index/{sample}.log",
    conda:
        "../envs/map_coverage_bin.yaml",
    benchmark: 
        "benchmarks/{sample}.indexCheckm.benchmark.txt",
    params:
        prefix = "../results/megahit_assm/{sample}_assm/{sample}",
        algorithm = "bwtsw",
    wrapper: 
        "0.80.2/bio/bwa/index"


#-------------------------------------------------------------------#
ruleorder: IndexMegahitContigs > MapReadsToContigs


#-------------------------------------------------------------------#
rule MapReadsToContigs: 
    """
        Map original reads to contigs
    """
    output: 
        sorted_bam = "../results/binning/{sample}/contig_index/{sample}_sorted.bam",
    input: 
        reads = megahit_input_1 and megahit_input_2,
        idx = "../results/megahit_assm/{sample}_assm/{sample}.sa",
    log: 
        "logs/bwa_mem/{sample}.log",
    threads: 
        25
    benchmark: 
        "benchmarks/{sample}.readMapping.benchmark.txt",
    params: 
        index = "../results/megahit_assm/{sample}_assm/{sample}",
        sorting = "samtools", 
        sort_order = 'coordinate',
    wrapper:
       "v1.17.3/bio/bwa/mem"


#-------------------------------------------------------------------#
rule MappingStats: 
    """
        Get mapping statistics using samtools
    """
    output: 
        map_stat = "../results/binning/{sample}/map_stat/{sample}_sorted.bam.stat",
    input: 
        bam_idx = "../results/binning/{sample}/contig_index/{sample}_sorted.bam",
    conda:
        "envs/map_coverage_bin.yaml",
    benchmark: 
        "benchmarks/{sample}.contigstat.benchmark.txt",
    shell: 
        r"""
            samtools flagstat {input.bam_idx} > {output.map_stat}
         """


#-------------------------------------------------------------------#
rule ConvertBamToDepth: 
    """
        Generate depth file of mapped reads 
    """
    output: 
        contig_depth = "../results/binning/{sample}/contig_depth/{sample}_sorted.depth",
    input: 
        bam_idx = "../results/binning/{sample}/contig_index/{sample}_sorted.bam",
    conda:
        "../envs/map_coverage_bin.yaml",
    benchmark: 
        "benchmarks/{sample}.contig_depth.benchmark.txt",
    shell: 
        r"""
            jgi_summarize_bam_contig_depths \
            --outputDepth {output.contig_depth} \
            {input.bam_idx}

         """


#-------------------------------------------------------------------#
rule MetabatBin: 
    """
        Bin assembled contigs using Metabat2
    """
    output: 
        bins = "../results/binning/metabat_out/{sample}_bins/{sample}",
    input: 
        contigs = "../results/megahit_assm/{sample}_assm/final.contigs.fa",
        contig_depth = "../results/binning/{sample}/contig_depth/{sample}_sorted.depth", 
    conda:
        "../envs/map_coverage_bin.yaml",
    benchmark: 
        "benchmarks/{sample}.metabat.benchmark.txt",
    params: 
        outdir = "../results/binning/metabat_out/{sample}_bins/bin/",
        min_contig_len = config["MetagenomeBinning"]["MinimumContigLength"],
        threads = config["MetagenomeBinning"]["Threads"],
    shell: 
        r"""
            metabat2 -i {input.contigs} \
            -a {input.contig_depth} \
            -o {output.bins} \
            -t {params.threads} \
            -m {params.min_contig_len} \
            --maxEdges 500

            touch {output}

         """

