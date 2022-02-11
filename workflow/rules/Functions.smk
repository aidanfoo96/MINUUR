"""
    Functions file. Series of statments to run the whole pipeline
    Determined by "Activate" status in configuration file
"""

def RunPipeline(wildcards):
    
    final_input = []
    
    if config["RemoveHostFromFastqGz"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/aligned_bam/{sample}_sorted.bam",
                ],
                sample = samples
            )
        )
    
    if config["QC"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/qc/fastqc/{sample}_{num}_fastqc.zip",
                    "../results/qc/trimmed_fastq/{sample}_trimmed_1.fastq",
                    "../results/qc/trimmed_fastq/{sample}_trimmed_2.fastq",
                    "../results/qc/fastqc_trimmed/{sample}_{num}.html", 

                ],
                sample = samples,
                num = [1, 2],
            )
        )

    if config["ProcessBam"]["Activate"]: 
        if config["ProcessBam"]["FromFastq"]:
            final_input.extend(
                expand(
                    [
                        "../results/bam_stats_ffq/{sample}_stats.txt",
                        "../results/unmapped_fastq_ffq/{sample}_unmapped_1.fastq",
                        "../results/alignment_stats_ffq/{sample}_align_stats.txt",
                        "../results/alignment_stats_ffq/concatenated_alignment_statistics.txt",
                        "../results/unmapped_bam_ffq/{sample}_unmapped.bam",
                    ],
                    sample = samples
                )
            )

        else:
            final_input.extend(
                expand(
                    [
                        "../results/bam_stats/{sample}_stats.txt",
                        "../results/alignment_stats/{sample}_align_stats.txt",
                        "../results/alignment_stats/concatenated_alignment_statistics.txt",
                        "../results/unmapped_bam/{sample}_unmapped.bam",
                        "../results/unmapped_fastq/{sample}_unmapped_1.fastq",
                    ],
                    sample = samples
                )
            )
        
    
    if config["KrakenClassification"]["Activate"]:
        final_input.extend(
            expand(
                [
                    "../results/kraken_out/mpa_report/{sample}_report.txt",
                    "../results/kraken_out/report/{sample}_kraken_report",
                ],
                sample = samples
            )
        )
                
    if config["KrakenSummaries"]["Activate"]:
        final_input.extend(
            expand(
                [
                    "../results/kraken_out/mpa_out/{sample}_mpa_conv_report.txt",
                    "../results/kraken_results/tables/kingdom_table_tidy.txt",
                    "../results/kraken_results/tables/genus_table_tidy.txt",
                    "../results/kraken_results/tables/species_table_tidy.txt",
                    "../results/kraken_results/tables/classified_reads_table.txt",
                    "../results/kraken_results/plots/classified_reads_plot.pdf",
                    "../results/alignment_stats_ffq/tables/aligned_reads_wide.txt",
                    "../results/alignment_stats_ffq/tables/aligned_reads_long.txt",
                    "../results/kraken_results/plots/classified_proportions.pdf",
                    "../results/alignment_stats_ffq/plots/alignment_stat_plot.pdf",
                    "../results/kraken_results/plots/genus_heatmap.pdf", 
                    "../results/kraken_results/classified_summary/{sample}_classified_summary.txt", 
                    "../results/kraken_results/concatenated_kraken_summary.txt",
                    "../results/kraken_results/plots/genus_spatial_plot.pdf", 
                    "../results/kraken_results/plots/species_spatial_plot.pdf",
                ],
                sample = samples
            )
        )

    if config["ExtractKrakenTaxa"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/kraken_taxon_extract/{sample}_extract_1.fastq",
                    "../results/kraken_taxon_extract/{sample}_extract_2.fastq",
                ],
                sample = samples
            )
        )

    if config["BrackenReestimation"]["Activate"]:
        final_input.extend(
            expand(
                [
                    "../results/bracken_reestimation/bracken_out/{sample}_bracken.txt", 
                    "../results/bracken_reestimation/concat_bracken_out/concatenated_bracken_report.txt",
                    "../results/bracken_reestimation/plots/stratified_species_heatmaps.pdf",
                ],
                sample = samples
            )
        )

    if config["MetaphlanClassification"]["Activate"]:
        final_input.extend(
            expand(
                [
                    "../results/metaphlan_out/taxa_profile/{sample}_taxa_prof.txt",
                    "../results/metaphlan_out/bowtie2_aln/{sample}.bowtie2.bz2",
                ],
                sample = samples
            )
        )

    if config["MetaphlanClassification"]["CleanMetaphlanReport"]: 
        final_input.extend(
            expand(
                [
                    "../results/metaphlan_out/taxa_profile_clean/{sample}_taxa_prof_clean.txt",
                    "../results/metaphlan_out/taxa_profile_clean/merged_tbl/merged_taxa_prof_clean.txt",
                    "../results/metaphlan_out/clean_summaries/kingdom_table_tidy.txt",
                    "../results/metaphlan_out/clean_summaries/genus_table_tidy.txt",
                    "../results/metaphlan_out/clean_summaries/species_table_tidy.txt"
                ],
                sample = samples
            )
        )

    if config["HumannAnalysis"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/concatenated_fastq/{sample}_unmapped_conc.fastq",
                    "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_genefamilies.tsv",
                    "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_pathabundance.tsv",
                    "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_pathcoverage.tsv",
                    "../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_genefamilies.tsv", 
                    "../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_pathabundance.tsv", 
                    "../results/humann_out/sorted_abundance_profiles/{sample}_unmapped_conc_pathcoverage.tsv",
                    "../results/humann_out/concatenated_humann_files/all_samples_genesfamilies_humann3.tsv",
                    "../results/humann_out/concatenated_humann_files/all_samples_pathabundance_humann3.tsv", 
                    "../results/humann_out/concatenated_humann_files/all_samples_pathcoverage_humann3.tsv",
                    "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_genefamilies_humann3_relab.tsv",
                    "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_pathabundance_humann3_relab.tsv",
                    "../results/humann_out/concatenated_humann_files/relative_abund/all_samples_pathcoverage_humann3_relab.tsv",
                    "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_genefamilies_humann3_relab_stratified.tsv",
                    "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_genefamilies_humann3_relab_unstratified.tsv",
                    "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathabundance_humann3_relab_stratified.tsv",
                    "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathabundance_humann3_relab_unstratified.tsv",
                    "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathcoverage_humann3_relab_stratified.tsv",
                    "../results/humann_out/concatenated_humann_files/relab_strat/all_samples_pathcoverage_humann3_relab_unstratified.tsv",
                    "../results/humann_out/summarised_bowtie2_stats/{sample}_bowtie2_alignment_summarised_gene_number.tsv",
                    "../results/humann_out/summarised_bowtie2_stats/concatenated_bowtie2_alignment_summarised_gene_number.tsv",
                    "../results/humann_out/{sample}_humann3_profile/{sample}_unmapped_conc_humann_temp/{sample}_unmapped_conc_bowtie2_aligned.tsv",
                ],
                sample = samples
            )
        )

    if config["HumannAnalysis"]["Activate"]:
        final_input.extend(
            expand(
                [
                    "../results/humann_out/concatenated_humann_files/renamed_abund_tbl/all_samples_genesfamilies_humann3_relab_stratified_uniref90.tsv",
                ]
                
            )
        )

    if config["GetBiologicalProcess"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/humann_out/extracted_function/{role}_extracted_function_bar.pdf",
                    "../results/humann_out/extracted_function/{role}_extracted_function_heatmap.pdf",
                ],
                role = config["GetBiologicalProcess"]["Process"]
            )
        )
    
    if config["MetagenomeAssm"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/megahit_assm/{sample}_assm/final.contigs.fa",
                    "../results/quast_out/{sample}_quast_output/report.txt",
                    "../results/quast_out/concat_transposed_report.tsv",
                    


                ],
                sample = samples

            )
        )


    if config["MetagenomeBinning"]["Activate"]: 
        final_input.extend(
            expand(
                [

                    "../results/binning/{sample}/contig_index/{sample}_sorted.bam",
                    "../results/binning/{sample}/contig_depth/{sample}_sorted.depth",
                    "../results/binning/metabat_out/{sample}_bins/{sample}",
                    "../results/megahit_assm/{sample}_assm/{sample}.amb",
                    "../results/megahit_assm/{sample}_assm/{sample}.ann",
                    "../results/megahit_assm/{sample}_assm/{sample}.bwt",
                    "../results/megahit_assm/{sample}_assm/{sample}.pac",
                    "../results/megahit_assm/{sample}_assm/{sample}.sa",
                    
                ],
                
                sample = samples

            )
        )   

    if config["CheckmBinQA"]["Activate"]: 
        final_input.extend(
            expand(
                [
                    "../results/binning/checkm_out/{sample}_checkm/lineage.ms",
                    "../results/binning/checkm_out/QA_out/{sample}_checkm_out.tsv",
                    "../results/binning/checkm_out/concatenated_bin_stats.tsv",
                    "../results/binning/plots/CompletenessVsContam.pdf",
                    "../results/binning/plots/NumContigsVsCompleteness.pdf",
                    "../results/binning/plots/CompletenessVsContam_Per_Sample.pdf",
                    "../results/binning/plots/BarChartCompletenessContamination.pdf", 

                ],
                sample = samples
            )
        )

    return(final_input)
