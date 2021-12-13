################ Quality Control and Host Alignment ##############
    # fastqc: (Andrews et al., 2010)
    # cutadapt: (Martin et al., 2011)
    # hisat2: (Kim et al., 2019)
##################################################################

def get_input(wildcards): 
    """
        Function to determine if input is raw fastq or BAM
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
    

rule do_fastqc: 
    """
        Perform FASTQC on paired fastq files
    """
    output: 
        html = "../results/qc/fastqc/{sample}_{num}.html", 
        zip = "../results/qc/fastqc/{sample}_{num}_fastqc.zip"
    input: 
        get_input
    log: 
        "logs/fastqc/{sample}_{num}.log", 
    threads: 
        2
    params: 
        outdir="--outdir ../results/qc/fastqc",
    wrapper: 
        "v0.80.1/bio/fastqc"


rule trim_sequences: 
    """
        Trim sequences fastq files 
        Trimming parameters specified in contig file under "CutadaptParams" 
    """
    output: 
        fastq1 = "../results/qc/trimmed_fastq/{sample}_trimmed_1.fastq",
        fastq2 = "../results/qc/trimmed_fastq/{sample}_trimmed_2.fastq", 
        qc="trimmed/{sample}.qc.txt",
    input:
        get_input
    params: 
        extra = config["QC"]["CutadaptParams"]
    threads: 
        4
    wrapper:
        "0.77.0/bio/cutadapt/pe"

rule do_fastqc_trimmed: 
    """
        Perform fastqc again on trimmed fastq files
    """
    output: 
        html = "../results/qc/fastqc_trimmed/{sample}_{num}.html", 
        zip = "../results/qc/fastqc_trimmed/{sample}_{num}_fastqc.zip"
    input: 
        fastq1 = "../results/qc/trimmed_fastq/{sample}_trimmed_1.fastq",
        fastq2 = "../results/qc/trimmed_fastq/{sample}_trimmed_2.fastq", 
    log: 
        "logs/fastqc_trimmed/{sample}_{num}.log", 
    threads: 
        2
    params: 
        outdir="--outdir ../results/qc/fastqc_trimmed",
    wrapper: 
        "v0.80.1/bio/fastqc"


rule align_fastq:
    """ 
        Align timmed fastq files against user defined indexed reference genome
        Generates coordinate sorted bam file
    """
    output: 
        bam = "../results/aligned_bam/{sample}_sorted.bam",
    input:
        read1 = "../results/qc/trimmed_fastq/{sample}_trimmed_1.fastq",
        read2 = "../results/qc/trimmed_fastq/{sample}_trimmed_2.fastq", 
    conda: 
        "envs/qc_env.yaml"
    params: 
        db = config['RemoveHostFromFastqGz']["ContaminantIndex"],
        threads = config['RemoveHostFromFastqGz']['Threads'],
    shell: 
        r"""
            hisat2 -x {params.db} \
            -1 {input.read1} -2 {input.read2} \
            -t --threads {params.threads} | 
            samtools view -bS - | samtools sort - \
            -o {output.bam}
         """
