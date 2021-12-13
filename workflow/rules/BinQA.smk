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
    params: 
        outdir = "../results/binning/checkm_out/{sample}_checkm/",
        indir = "../results/binning/metabat_out/{sample}_bins/",
        threads = config["CheckmBinQA"]["Threads"],
    shell: 
        r"""
            checkm lineage_wf -t {params.threads} \
            --pplacer_threads {params.threads} \
            -x .fa {params.indir} {params.outdir}    
         """

rule QA_checkm: 
    """
        Run qa from checkm to generate .tsv summary file
    """
    output:
        "../results/binning/checkm_out/QA_out/{sample}_checkm_out.tsv",
    input: 
        lineage = "../results/binning/checkm_out/{sample}_checkm/lineage.ms",
    params: 
        threads = config["CheckmBinQA"]["Threads"],
        out = "../results/binning/{sample}/checkm_out/lineage.ms",
        direc = "../results/binning/checkm_out/{sample}_checkm/",
    shell: 
        r"""
            checkm qa -t {params.threads} -f {output} --tab_table -o 2 \
            {input.lineage} {params.direc}
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
    input: 
        checkm_conc_result = "../results/binning/checkm_out/concatenated_bin_stats.tsv",
    script:
        "../scripts/plot_checkm_QA.R" 