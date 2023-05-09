################ Bin Quality Assurance ##############
    # checkm: (Parks et al., 2015)
#####################################################

rule lineage_checkm: 
    """
        Run lineage_wf from checkm to assess quality of metabat bins
    """
    output: 
        "../results/binning/checkm_out/{sample}_checkm/lineage.ms",
    input: 
        direc = "../results/binning/metabat_out/{sample}_bins/{sample}",
    log: 
        "logs/BinQA/{sample}.log",
    conda:
        "../envs/map_coverage_bin.yaml",
    benchmark: 
        "benchmarks/{sample}.checkMLineage.benchmark.txt",
    params: 
        outdir = "../results/binning/checkm_out/{sample}_checkm/",
        indir = "../results/binning/metabat_out/{sample}_bins/",
        threads = config["CheckmBinQA"]["Threads"],
    shell: 
        r"""
            checkm lineage_wf -t {params.threads} \
            --pplacer_threads {params.threads} \
            -x .fa {params.indir} {params.outdir} 2> {log}  
         """

rule QA_checkm: 
    """
        Run qa from checkm to generate .tsv summary file
    """
    output:
        "../results/binning/checkm_out/QA_out/{sample}_checkm_out.tsv",
    input: 
        lineage = "../results/binning/checkm_out/{sample}_checkm/lineage.ms",
    log: 
        "logs/checkmQA/{sample}.log"
    params: 
        threads = config["CheckmBinQA"]["Threads"],
        out = "../results/binning/{sample}/checkm_out/lineage.ms",
        direc = "../results/binning/checkm_out/{sample}_checkm/",
    benchmark: 
        "benchmarks/{sample}.checkMQA.benchmark.txt",
    conda:
        "../envs/map_coverage_bin.yaml",
    shell: 
        r"""
            checkm qa -t {params.threads} -f {output} --tab_table -o 2 \
            {input.lineage} {params.direc} 2> {log}
         """

rule concatenate_checkm_out: 
    """
        Concatenate all qa .tsv files from each checkm run per sample
    """
    output: 
        txt = "../results/binning/checkm_out/concatenated_bin_stats.tsv",
    input: 
        sample = expand("../results/binning/checkm_out/QA_out/{sample}_checkm_out.tsv", sample = samples),
    params: 
        filename = "FILENAME", 
    shell: 
        r"""
            awk '{{print $0 "\t" {params.filename}}}' {input.sample} > {output.txt}
         """

rule plot_QA_stats: 
    """
        Use concatenated .tsv file to generate summary plots 
        1. Completeness vs Contamination Bubble plot
        2. Number of Contigs vs Completeness Plot
    """
    output: 
        CompletenessVsContam = "../results/binning/plots/CompletenessVsContam.pdf",
        NumContigsVsCompleteness = "../results/binning/plots/NumContigsVsCompleteness.pdf",
        CompletenessVsContam_Per_Sample = "../results/binning/plots/CompletenessVsContam_Per_Sample.pdf",
        BarChartCompletenessContamination = "../results/binning/plots/BarChartCompletenessContamination.pdf", 
    input: 
        checkm_conc_result = "../results/binning/checkm_out/concatenated_bin_stats.tsv",
    conda: 
        "../envs/r_and_plotting_env.yaml",
    script:
        "../scripts/plot_checkm_QA.R" 

rule remove_NoExtensionMAGs: 
    """
        Removes MAGs that do not end in .fa - 0Kb
    """
    input: 
        data = expand("../results/binning/metabat_out/{sample}_bins/", sample = samples),
    shell: 
        r'''
            find {input.data} -not -name '*.fa' -type f -delete
        '''

rule run_BUSCO: 
    """
        Run BUSCO
        Checks is certain MAGs are eukaryotic
    """
    output:
        dataset_dir=directory("../results/binnning/metabat_out/busco/busco_downloads"),
    input: 
        direc = expand("../results/binning/metabat_out/{sample}_bins/", sample=samples),
    params:
        mode="genome",
        lineage="metazoa_odb10",
    threads: 8
    wrapper:
        "v1.28.0/bio/busco"