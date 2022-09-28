################ Quality Control and Host Alignment ##############
    # fastqc: (Andrews et al., 2010)
    # cutadapt: (Martin et al., 2011)
    # bowtie2: (Author?)
##################################################################

from snakemake.utils import validate

#-------------------------------------------------------------------#
def GetInput(wildcards): 
    """
        Function to determine if input is raw fastq or BAM
    """
    if config["input_config"]["automatic"]:
        units = pd.read_csv(config["samples"], sep = "\t")
        units = (
            units.assign(fastq1Path=f"../resources/" + units["sampleID"] + "_1.fastq.gz")
            .assign(fastq2Path=f"../resources/" + units["sampleID"] + "_2.fastq.gz")
            .set_index("sampleID")
        )
    else:
        assert os.path.isfile(
            config["input_config"]["manual"]
        ), f"config['input_config']['manual'] doesn't exist. Please create one or use the 'automatic option' with the correct naming scheme specified in the WIKI page"
        units = pd.read_csv(config["input_config"]["manual"], sep = "\t", index_col="sampleID")

    u = units.loc[wildcards.sample, ["fastq1Path", "fastq2Path"]].dropna()

    return [f"{u.fastq1Path}", f"{u.fastq2Path}"]

#-------------------------------------------------------------------#
rule CheckInputs:
    """
    Check all inputs are present 
    """
    input:
        ref = config["RemoveHostFromFastqGz"]["ContaminantIndex"],
    output:
        touch("../results/.input.check"),
    params:
        sample_list = config["samples"],
        sample_table = config["input_config"]["manual"],
        automatic_input = config["input_config"]["automatic"]
    log:
        "logs/CheckInputs.log",
    priority: 50
    script: 
        "../scripts/check_inputs.py"


validate(sample_table, schema = "../schemas/samples.schema.yaml")

#-------------------------------------------------------------------#
rule FastQC: 
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
        10
    benchmark: 
        "benchmarks/{sample}_{num}.FastQC.benchmark.txt",
    params: 
        outdir="--outdir ../results/qc/fastqc",
    wrapper: 
        "v0.80.1/bio/fastqc"

#-------------------------------------------------------------------#
rule TrimFastq: 
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
        extra = config["QC"]["CutadaptParams"],
    log: 
        "logs/cutadapt/{sample}.log", 
    benchmark: 
        "benchmarks/{sample}.trimming.benchmark.txt",
    wrapper:
        "0.77.0/bio/cutadapt/pe"


#-------------------------------------------------------------------#
rule FastQCTrimmed: 
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
    benchmark: 
        "benchmarks/{sample}_{num}.FASTQC_Trimmed.benchmark.txt",
    threads: 
        10
    params: 
        outdir="--outdir ../results/qc/fastqc_trimmed",
    wrapper: 
        "v0.80.1/bio/fastqc"


#-------------------------------------------------------------------#
rule AlignFastq:
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
        "../envs/qc_env.yaml",
    log: 
        "logs/bowtie2_align/{sample}.log", 
    benchmark: 
        "benchmarks/{sample}.alignment.benchmark.txt",
    params: 
        db = config['RemoveHostFromFastqGz']["ContaminantIndex"],
        threads = config['RemoveHostFromFastqGz']['Threads'],
        sensitivity = config['RemoveHostFromFastqGz']['AlignmentSensitivty']
    shell: 
        r"""
            bowtie2 -x {params.db} \
            -1 {input.read1} -2 {input.read2} \
            -t --threads {params.threads} {params.sensitivity} | 
            samtools view -bS - | samtools sort - \
            -o {output.bam} 2> {log}
         """