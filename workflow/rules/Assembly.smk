################ Metagenome Assembly and Statistics ################
    # MegaHit: (Li et al., 2015) 
    # QUAST: (Gurevich et al., 2013) 
####################################################################

# Functions to return fastq files depending on configuration 
def megahit_input_1(wildcards): 
    if config["MetagenomeAssm"]["UseKrakenExtracted"]:             
        return("../results/kraken_taxon_extract/{sample}_extract_1.fastq")
    elif config["input_type"]["fastq"]:
        return("../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq")
    elif config["input_type"]["bam"]: 
        return("../results/unmapped_fastq/{sample}_unmapped_1.fastq")
    

def megahit_input_2(wildcards): 
    if config["MetagenomeAssm"]["UseKrakenExtracted"]: 
        return("../results/kraken_taxon_extract/{sample}_extract_2.fastq")
    elif config["input_type"]["fastq"]:
        return("../results/unmapped_fastq_ffq/{sample}_unmapped_2.fastq")
    elif config["input_type"]["bam"]: 
        return("../results/unmapped_fastq/{sample}_unmapped_2.fastq")


rule megahit_assembly: 
    """
        Perform metagenome assembly with MegaHit
        Inputs = unmapped fastq OR kraken extracted reads
    """
    output: 
        final_contigs = "../results/megahit_assm/{sample}_assm/final.contigs.fa"
    input: 
        read1 = megahit_input_1,
        read2 = megahit_input_2
    conda: 
        "../envs/metagenome_assm_env.yaml",
    log:
        "logs/megahit_assm/{sample}.log",
    params: 
        outdir = "../results/megahit_assm/{sample}_assm/tmp",
        threads = config["MetagenomeAssm"]["Threads"],
        memory = config["MetagenomeAssm"]["Memory"]
    shell: 
        r"""

            set +o pipefail 
            
            megahit \
            -1 {input.read1} \
            -2 {input.read2} \
            -o {params.outdir} \
            -t {params.threads} -m {params.memory} \
            --k-step 10 2> {log}

            mv ../results/megahit_assm/{wildcards.sample}_assm/tmp/* \
            ../results/megahit_assm/{wildcards.sample}_assm/
         """        


rule quast: 
    """
        Get contig statistics using QUAST
    """
    output: 
        quast_out = "../results/quast_out/{sample}_quast_output/report.txt",
    input: 
        megahit_out = "../results/megahit_assm/{sample}_assm/final.contigs.fa",
    conda: 
        "../envs/metagenome_assm_env.yaml",
    log:
        "logs/quast/{sample}.log",
    params: 
        outdir = "../results/quast_out/{sample}_quast_output/",
        threads = config["MetagenomeAssm"]["Threads"],
    shell: 
        r"""
            quast.py -o {params.outdir} \
            -t {params.threads} \
            --labels megahit \
            {input.megahit_out} 2> {log}

         """


rule concatenate_quast_report: 
    """
        Combine QUAST reports into one file for parsing
    """
    output:    
        txt =  "../results/quast_out/concat_transposed_report.tsv"
    input: 
        sample = expand("../results/quast_out/{sample}_quast_output/transposed_report.tsv", sample = samples), 
    params: 
        filename = "FILENAME", 
    shell: 
        r"""
            awk '{{print $0 "\t" {params.filename}}}' {input.sample} > {output.txt}
         """
