#!/usr/bin/env Rscript
library(tidyverse)

# Load Appropriate Columns 
columns <- c("bin", "marker_lineage", "num_genomes", "num_markers", "num_marker_sets", "completeness", "contamination",
             "strain_heterogen", "genome_size_bp", "num_ambiguous_bases", "num_scaffolds", "num_contigs", "N50_scaffold", "N50_contig",
             "mean_scaffold_length", "mean_contig_length", "longest_scaffold", "longest_contig", 
             "GC", "GC_std", "coding_density", "translation_table", "num_predicted_genes", "0", "1", "2", "3", "4", "5+", "sample")

# Read Data, skip first row and remove the redundant header columns from previous step
bin_stats <- read_tsv(snakemake@input[["checkm_conc_result"]], col_names = columns, skip = 1) %>%
  filter(bin != "Bin Id")

# Conver to Numeric Columns
bin_stats$completeness <- as.numeric(bin_stats$completeness)
bin_stats$contamination <- as.numeric(bin_stats$contamination)
bin_stats$genome_size_bp <- as.numeric(bin_stats$genome_size_bp)
bin_stats$num_contigs <- as.numeric(bin_stats$num_contigs)
bin_stats$mean_contig_length <- as.numeric(bin_stats$mean_contig_length)
bin_stats$N50_contig <- as.numeric(bin_stats$N50_contig)

# Get breaks for plotting 
bin_stats$col <- cut(bin_stats$completeness, 
                        breaks = c(-Inf, 90, Inf), 
                        labels = c("<90", ">=90"))
  
bin_stats$contig_cut <- cut(bin_stats$num_contigs, 
                     breaks = c(-Inf, 250, Inf), 
                     labels = c("<=250", ">250"))

CompletenessVsContam <- bin_stats %>%
  ggplot() + 
  aes(x = completeness, y = contamination, size = genome_size_bp, col = col) %>%
  geom_point(alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 100 )) +
  theme_minimal() + 
  scale_color_manual(values = c("blue", "orange"))

pdf(snakemake@output[["CompletenessVsContam"]])
print(CompletenessVsContam)
dev.off()


NumContigsVsCompleteness <- bin_stats %>%
  ggplot() + 
  aes(x = num_contigs, y = completeness, size = mean_contig_length, col = contig_cut) + 
  geom_point(alpha = 0.5) + 
  theme_minimal() + 
  scale_color_manual(values = c("blue", "orange"))

pdf(snakemake@output[["NumContigsVsCompleteness"]])
print(NumContigsVsCompleteness)
dev.off()

