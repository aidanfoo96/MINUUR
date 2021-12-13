################ Metagenome Co-Assembly ##############
    # IN PROGRESS NOT FINISHED 
#####################################################

def get_co_assm_input1(wildcards): 

    input_list = []

    if config["MetagenomeCoAssm"]["Activate"]: 
        if config["MetagenomeAssm"]["UseKrakenExtracted"]:
            input_list.extend(
                expand(
                    [
                        "../results/kraken_taxon_extract/{sample}_extract_1.fastq",
                    ],
                    sample = samples,
                )
            )
        else: 
            input_list.extend(
                expand(
                    [
                        "../results/unmapped_fastq/{sample}_unmapped_1.fastq",
                    ],
                    sample = samples, 
                )
            )
    
        return(input_list)

def get_co_assm_input2(wildcards): 

    input_list = []

    if config["MetagenomeAssm"]["MetagenomeCoAssm"]: 
        if config["MetagenomeAssm"]["UseKrakenExtracted"]:
            input_list.extend(
                expand(
                    [
                        "../results/kraken_taxon_extract/{sample}_extract_2.fastq",
                    ],
                    sample = samples,
                )
            )
        else: 
            input_list.extend(
                expand(
                    [
                        "../results/unmapped_fastq/{sample}_unmapped_2.fastq",
                    ],
                    sample = samples, 
                )
            )
    
        return(input_list)


rule co_assembly: 
    output:
        final_contigs = "../results/megahit_coassm/final.contigs.fa", 
    input: 
        read1 = get_co_assm_input1,
        read2 = get_co_assm_input2,
    params: 
        outdir = "../results/megahit_coassm/tmp", 
        threads = config["MetagenomeCoAssm"]["Threads"],
        memory = config["MetagenomeCoAssm"]["Memory"],
    shell: 
        r"""
            set +o pipefail 
            
            megahit \
            -1 {input.read1} \
            -2 {input.read2} \
            -o {params.outdir} \
            -t {params.threads} -m {params.memory} 

            mv ../results/megahit_coassm/tmp/* \
            ../results/megahit_coassm/
         """        

    