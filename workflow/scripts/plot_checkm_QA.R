#!/usr/bin/env Rscript
library(tidyverse)
library(MetBrewer)

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
                        breaks = c(-Inf, 80, Inf), 
                        labels = c("<80", ">=80"))
  
bin_stats$contig_cut <- cut(bin_stats$num_contigs, 
                     breaks = c(-Inf, 250, Inf), 
                     labels = c("<=250", ">250"))

# Plot Completeness Vs Contamination Per Sample
CompletenessVsContam <- bin_stats %>%
  ggplot() + 
  aes(x = completeness, y = contamination, size = genome_size_bp, col = col) %>%
  geom_point(alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 100)) +
  theme_minimal() + 
  scale_color_manual(values = c("blue", "orange"))

pdf(snakemake@output[["CompletenessVsContam"]])
print(CompletenessVsContam)
dev.off()

#### Plot Completeness Vs Contamination Per Sample
bin_stats_sample_split <- bin_stats %>%
  separate(bin, into = c("sample", "bin_num"), sep = "\\.")

compvscontam_samp <- bin_stats_sample_split %>%
  ggplot() + 
  aes(x = completeness, y = contamination, col = sample) %>%
  geom_point(alpha = 0.8, size = 4) + 
  scale_x_continuous(limits = c(0, 100 )) +
  geom_vline( xintercept = 80, linetype = "dotted", 
              color = "green", size = 1) + 
  geom_vline(xintercept = 50, linetype = "dotted", 
             color = "orange", size = 1) + 
  theme_bw(base_size = 18) + 
  scale_color_manual(values = met.brewer("Redon", 10, type = "discrete")) + 
  xlab("Completeness") + 
  ylab("Contamination") + 
  labs(colour = "Sample") + 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.text = element_text(size = 10))

pdf(snakemake@output[["CompletenessVsContam_Per_Sample"]])
print(compvscontam_samp)
dev.off()

#### Plot Bar Chart with Completeness vs Contamination Per Sample
bin_stats_bar <- bin_stats_sample_split %>%
  unite("sample_bin", sample:bin_num) %>%
  select(completeness, contamination, sample_bin) %>%
  pivot_longer(c(completeness, contamination), 
              names_to = "bin_rank", 
              values_to = "value") %>%
  ggplot() + 
  aes(reorder(x = sample_bin, desc(value)), y = value, fill = bin_rank) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw(base_size = 15) + 
  coord_flip() + 
  scale_fill_manual(values = met.brewer("Redon", 2, type = "discrete"), labels = c("Completeness", "Contamination")) + 
  geom_hline(yintercept = 80, linetype = "dotted", 
              color = "green", size = 1) + 
  geom_hline(yintercept = 50, linetype = "dotted", 
             color = "orange", size = 1) + 
  xlab("Sample") + 
  ylab("Contamination and Completion Score (%)") + 
  labs(fill = "CheckM Completion vs Contamination Score") + 
  theme(legend.position = "bottom", 
        legend.direction = "vertical")

pdf(snakemake@output[["BarChartCompletenessContamination"]])
print(bin_stats_bar)
dev.off()

#### Plot the number of contigs vs completeness
NumContigsVsCompleteness <- bin_stats %>%
  ggplot() + 
  aes(x = num_contigs, y = completeness, size = mean_contig_length, col = contig_cut) + 
  geom_point(alpha = 0.5) + 
  theme_minimal() + 
  scale_color_manual(values = c("blue", "orange"))

pdf(snakemake@output[["NumContigsVsCompleteness"]])
print(NumContigsVsCompleteness)
dev.off()

