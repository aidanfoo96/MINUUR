################ Functional Read Profiling #####################
    # Humann3: (Beghini et al., 2021)
###################################################################


def get_unmap_fastq_1(wildcards): 
    if config["input_type"]["fastq"]: 
        return("../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq")
    elif config["input_type"]["bam"]: 
        return("../results/unmapped_fastq/{sample}_unmapped_1.fastq")


def get_unmap_fastq_2(wildcards): 
    if config["input_type"]["fastq"]: 
        return("../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq")
    elif config["input_type"]["bam"]: 
        return("../results/unmapped_fastq/{sample}_unmapped_1.fastq")


rule concatenate_paired_reads: 
    """
        Concatenate unmapped fastq files 
    """
    output:
        concat_reads = "../results/concatenated_fastq/{sample}_unmapped_conc.fastq",
    input: 
        fq1 = get_unmap_fastq_1,
        fq2 = get_unmap_fastq_2,
    shell: 
        "cat {input.fq1} {input.fq2} > {output.concat_reads}"


rule generate_humann_profile:
    """
        Run Humann3 to generate humann functional profiles 
    """
    output: 
        gene_fam = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_genefamilies.tsv",
        path_abun = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_pathabundance.tsv",
        path_cov = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_pathcoverage.tsv",
        unmapped_humann_tmp = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_humann_temp/{sample}_unmapped_conc_bowtie2_aligned.tsv",
    input: 
        concat_reads = "../results/concatenated_fastq/{sample}_unmapped_conc.fastq",
        mtphlan_prof = "../results/metaphlan_out/taxa_profile/{sample}_taxa_prof.txt",
    conda: 
        "../envs/readfunction_env.yaml",
    params: 
        nt_db = config["HumannAnalysis"]["NTDatabase"],
        prot_db = config["HumannAnalysis"]["ProtDatabase"],
        num_threads = config["HumannAnalysis"]["Threads"],
        outdir = directory("../results/humann_out/{sample}_humann3_profile/"),
    log:
        "log/humann_profiling/{sample}.log",
    shell: 
        r"""
            humann3 --nucleotide-database {params.nt_db} \
            --protein-database {params.prot_db} \
            --threads {params.num_threads} \
            --input {input.concat_reads} \
            --output {params.outdir} \
            --taxonomic-profile {input.mtphlan_prof} 2> {log}
         """


rule sort_files:   
    """
        Sort huamnn outputs
    """
    output:
        moved_genes = "../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_genefamilies.tsv", 
        moved_abun = "../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_pathabundance.tsv", 
        moved_coverage = "../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_pathcoverage.tsv",
    input: 
        gene_fam = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_genefamilies.tsv",
        path_abun = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_pathabundance.tsv",
        path_cov = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_pathcoverage.tsv",
    conda: 
        "../envs/readfunction_env.yaml",
    shell: 
        r"""
            cp {input.gene_fam} {output.moved_genes}
            cp {input.path_abun} {output.moved_abun}
            cp {input.path_cov} {output.moved_coverage}
         """


rule join_humann_outputs: 
    """
        Join humann outputs from all samples specified in run
    """
    output: 
        joined_genes = "../results/humann_out/concatenated_humann_files/all_samples_genesfamilies_humann3.tsv",
        joined_abundances = "../results/humann_out/concatenated_humann_files/all_samples_pathabundance_humann3.tsv", 
        joined_coverages = "../results/humann_out/concatenated_humann_files/all_samples_pathcoverage_humann3.tsv",
    input: 
        moved_genes = expand("../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_genefamilies.tsv", sample=samples),
        moved_abun = expand("../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_pathabundance.tsv", sample = samples),
        moved_coverage = expand("../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_pathcoverage.tsv", sample=samples),
    params: 
        indir = "../results/humann_out/sorted_abundance_profiles/",
    conda: 
        "../envs/readfunction_env.yaml",
    shell: 
        r"""
            humann_join_tables -i {params.indir} -o {output.joined_genes} --file_name genefamilies
            humann_join_tables -i {params.indir} -o {output.joined_abundances} --file_name pathabundance
            humann_join_tables -i {params.indir} -o {output.joined_coverages} --file_name pathcoverage
         """


rule normalise_abundance:
    """
        Normalise humann read levels 
    """
    output: 
        relab_genes = "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_genefamilies_humann3_relab.tsv",
        relab_path = "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_pathabundance_humann3_relab.tsv",
        relab_cov = "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_pathcoverage_humann3_relab.tsv",
    input: 
        joined_genes = "../results/humann_out/concatenated_humann_files/all_samples_genesfamilies_humann3.tsv",
        joined_abundances = "../results/humann_out/concatenated_humann_files/all_samples_pathabundance_humann3.tsv", 
        joined_coverages = "../results/humann_out/concatenated_humann_files/all_samples_pathcoverage_humann3.tsv",
    conda: 
        "../envs/readfunction_env.yaml",
    
    shell: 
        r"""
            humann_renorm_table -i {input.joined_genes} -u relab -p -o {output.relab_genes}
            humann_renorm_table -i {input.joined_abundances} -u relab -p -o {output.relab_path}
            humann_renorm_table -i {input.joined_coverages} -u relab -p -o {output.relab_cov}
         """


rule split_output: 
    """
        Generate stratified and unstratified humann outputs
    """
    output:
        g_strat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_genefamilies_humann3_relab_stratified.tsv",
        g_unstrat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_genefamilies_humann3_relab_unstratified.tsv",
        p_strat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathabundance_humann3_relab_stratified.tsv",
        p_unstrat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathabundance_humann3_relab_unstratified.tsv",
        c_strat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathcoverage_humann3_relab_stratified.tsv",
        c_unstrat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathcoverage_humann3_relab_unstratified.tsv",
    input: 
        relab_genes = "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_genefamilies_humann3_relab.tsv",
        relab_path = "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_pathabundance_humann3_relab.tsv",
        relab_cov = "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_pathcoverage_humann3_relab.tsv",
    params: 
        outdir = directory("../results/humann_out/concatenated_humann_files/relab_strat/"), 
    conda: 
        "../envs/readfunction_env.yaml",
    shell: 
        r"""
            humann_split_stratified_table -i {input.relab_genes} -o {params.outdir}
            humann_split_stratified_table -i {input.relab_path} -o {params.outdir}
            humann_split_stratified_table -i {input.relab_cov} -o {params.outdir}
         """


rule rename_humann_genes: 
    """ 
        Rename the humann outputs to a user specified database 
    """
    output:
        renamed_genes = "../results/humann_out/concatenated_humann_files/renamed_abund_tbl/all_samples_genesfamilies_humann3_relab_stratified_uniref90.tsv",
    input: 
        g_strat = "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_genefamilies_humann3_relab_stratified.tsv",
    params: 
        rename = config["RenameHumannGeneNames"]["Rename"],
    conda: 
        "../envs/readfunction_env.yaml",
    shell: 
        r"""
            humann_rename_table --input {input.g_strat} \
            --output {output.renamed_genes} \
            -c {params.rename}
         """


rule plot_function: 
    """
        Generate user specified plots of desired function
    """
    output: 
        stacked_bar_plot = "../results/humann_out/extracted_function/{role}_extracted_function_bar.pdf",
        heatmap_plot = "../results/humann_out/extracted_function/{role}_extracted_function_heatmap.pdf",
    input: 
        renamed_table = "../results/humann_out/concatenated_humann_files/renamed_abund_tbl/all_samples_genesfamilies_humann3_relab_stratified_uniref90.tsv",
    params:
        role = config["GetBiologicalProcess"]["Process"],
    conda: 
        "../envs/r_and_plotting_env.yaml",
    script: 
        "../scripts/extract_humann_function.R"

rule get_bowtie2_alignment_statistics:
    """
        Get bowtie2 alignment sensitvity
    """
    output:
        summarised_genes = "../results/humann_out/summarised_bowtie2_stats/{sample}_bowtie2_alignment_summarised_gene_number.tsv",
    input: 
        bowtie2_humann_align = "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_humann_temp/{sample}_unmapped_conc_bowtie2_aligned.tsv", 
    conda: 
        "../envs/r_and_plotting_env.yaml",
    script: 
        "../scripts/extract_humann_bowtie2_percentidentity.R"

rule combine_bowtie2_humann_stat: 
    """
        Combine bowtie2 alignment_sensitvity plots
    """
    output: 
        combined_summarised_genes = "../results/humann_out/summarised_bowtie2_stats/concatenated_bowtie2_alignment_summarised_gene_number.tsv",
    input: 
        summarised_genes = expand("../results/humann_out/summarised_bowtie2_stats/{sample}_bowtie2_alignment_summarised_gene_number.tsv", sample=samples),
    params: 
        filename = "FILENAME",
    shell: 
        r"""
            awk '{{print $0 "\t" {params.filename}}}' {input.summarised_genes} > {output.combined_summarised_genes}
         """


## Add a rule to colour a specific taxa + function in plot! 
