#!/usr/bin/env Rscript
library(tidyverse)
library(MetBrewer)

#### Plot humann statistics
concatenated_alignment_stat <- read_tsv(snakemake@input[["humann_alignment_stats"]], 
                                        col_names = c("percent_iden", "gene_number", "sample", "path")) %>%
  filter(percent_iden != "percentage_identity")

genes <- read_tsv(snakemake@input[["gene_fam_norelab"]])

pathabund <- read_tsv(snakemake@input[["path_fam_norelab"]])

#### Plot ANI per gene ####
concatenated_alignment_stat_clean <- concatenated_alignment_stat %>%
  separate(sample, into = c("junk", "sampleID"), sep = "profile/") %>%
  select(sampleID, percent_iden, gene_number)

concatenated_alignment_stat_clean$percent_iden <- as.numeric(concatenated_alignment_stat_clean$percent_iden)
concatenated_alignment_stat_clean$gene_number <- as.numeric(concatenated_alignment_stat_clean$gene_number)

plot_gene_percent_number <- function(data){
  plot <- data %>%
    ggplot() + 
    aes(x = sampleID, y = percent_iden, size = gene_number) + 
    geom_jitter(alpha = 0.3, width = 0.3) + 
    scale_size(range = c(.1, 15)) + 
    ylim(0, 100) + 
    theme_bw(base_size = 15) + 
    ylab("Percentage Identity (%)") + 
    xlab("Sample") + 
    theme(legend.position = "top") + 
    guides(size = guide_legend(title = "Number of Identified Genes"))
  return(plot)
}

#### Wrangle Gene and Pathway Data ####
genes_clean <- genes %>%
  separate(`# Gene Family`, c("gene", "genus"), sep = "g__") %>%
  separate(genus, c("genus", "species"), sep = "s__") %>%
  separate(gene, c("gene", "classification"), sep = "\\|") %>%
  filter(classification != "unclassified") %>%
  select(!classification)

genes_clean_long <- genes_clean %>%
  pivot_longer(cols = !c(gene, genus, species), 
               names_to = "sample", 
               values_to = "reads_RPK") %>%
  filter(reads_RPK > 0) 

pathabund_clean <- pathabund %>%
  separate(`# Pathway`, c("pathway", "genus"), sep = "g__") %>%
  separate(genus, c("genus", "species"), sep = "s__") %>%
  separate(pathway, c("label", "description"), sep = ":") %>%
  filter(label != "UNINTEGRATED|unclassified") %>%
  filter(label != "UNINTEGRATED|") %>%
  filter(label != "UNMAPPED") %>%
  filter(label != "UNINTEGRATED") 

pathabund_clean_long <- pathabund_clean %>%
  pivot_longer(cols = !c(label, description, species, genus), 
               names_to = "sample", 
               values_to = "path_abund") %>%
  filter(path_abund > 0) %>%
  filter(genus != "NA") %>%
  filter(description  != "NA")

#### Plot Gene and Pathway Data ####
plot_abundance_per_sample <- function(humann_genes_long, humann_paths_long, ylabel){
  
  genes_per_sample <- humann_genes_long %>%
    group_by(sample) %>%
    summarise(identified_genes = n()) %>%
    separate(sample, c("sample", "junk2"), sep = "_unmapped") 
  
  
  pathways_per_sample <- humann_paths_long %>%
    group_by(sample) %>%
    summarise(identified_pathways = n()) %>%
    separate(sample, c("sample", "junk2"), sep = "_unmapped") 
  
  
  joined_table <- genes_per_sample %>%
    left_join(pathways_per_sample, by = "sample") %>%
    pivot_longer(cols = c(identified_genes, identified_pathways), 
                 names_to = "identification", 
                 values_to = "number")
  
  
  plot <- joined_table %>%
    separate(sample, c("sample", "junk2"), sep = "_unmapped") %>%
    ggplot() + 
    aes(reorder(x = sample, number), y = number, fill = identification) + 
    geom_bar(stat = "identity", width = 0.6, position = "dodge") + 
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1), 
          legend.title = element_blank(), 
          legend.position = "top") + 
    scale_fill_manual(values = met.brewer("Redon", 2), labels = c("Gene", "Pathway")) +
    ylab(ylabel) + 
    xlab("Sample") 
  
  return(plot)
}


plot_genes_per_species <- function(humann_genes_long, humann_paths_long, ylabel){
  
  genes_per_sample <- humann_genes_long %>%
    group_by(species) %>%
    summarise(identified_genes = n()) 
  
  pathways_per_sample <- humann_paths_long %>%
    group_by(species) %>%
    summarise(identified_pathways = n()) 
  
  joined_table <- genes_per_sample %>%
    left_join(pathways_per_sample, by = "species") %>%
    pivot_longer(cols = c(identified_genes, identified_pathways), 
                 names_to = "identification", 
                 values_to = "number") %>%
    filter(species != "NA")
  
  plot <- joined_table %>%
  ggplot() +
  aes(reorder(x = species, number), y = number, fill = identification) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") + 
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top") + 
  scale_fill_manual(values = met.brewer("Redon", 2), labels = c("Gene", "Pathway")) +
  ylab(ylabel) + 
  xlab("Species") + 
  coord_flip()
  
  return(plot)
  
}

gene_percent_iden <- plot_gene_percent_number(concatenated_alignment_stat_clean)
hits_per_sample <- plot_abundance_per_sample(genes_clean_long, pathabund_clean_long, "Summed Hits from All Samples")
hits_per_species <- plot_genes_per_species(genes_clean_long,pathabund_clean_long, "Summed Hits from All Samples")

ggsave(snakemake@output[["gene_hits"]], gene_percent_iden)
ggsave(snakemake@output[["hits_per_sample"]], hits_per_sample)
ggsave(snakemake@output[["hits_per_species"]], hits_per_species)