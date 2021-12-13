#!/usr/bin/env Rscript
# Load Packages 
library(tidyverse)
# Get colour pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Set Read Filter
filter <- as.character(snakemake@params[["bracken_threshold"]])
filter <- as.numeric(filter)

#### Load Data #####
conc_data <- read_tsv(snakemake@input[["conc_input"]], col_names = c("species",
                                                                       "taxa_id", 
                                                                       "taxa_lvl", 
                                                                       "kraken_assigned_reads", 
                                                                       "added_reads", 
                                                                       "new_est_reads", 
                                                                       "fraction_total", 
                                                                       "file_path")) %>%
  filter(species != "name")

conc_data_clean <- conc_data %>%
  separate(file_path, c("file", "sample"), sep = "out/") %>%
  separate(sample, c("sample", "junk"), sep = "_bracken") %>%
  select(sample, species, taxa_id, taxa_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total) %>%
  mutate_at(c("taxa_id", "kraken_assigned_reads", "added_reads", "new_est_reads", "fraction_total"), as.numeric)
  
  
#### Added Reads Across Samples ####

sample_summary <- conc_data_clean %>%
  group_by(sample) %>%
  summarise(kraken_reads = sum(kraken_assigned_reads), 
            added_reads = sum(added_reads), 
            total_new_reads = sum(new_est_reads))

sample_summary_long <- sample_summary %>%
  pivot_longer(cols = !sample, 
               names_to = "classification_type", 
               values_to = "num_reads") 

sample_summary_plot <- sample_summary %>%
  pivot_longer(cols = !sample, 
               names_to = "classification_type", 
               values_to = "num_reads") %>%
  filter(classification_type != "total_new_reads") %>%
  ggplot() + 
  aes(reorder(x = sample, desc(num_reads)), y = num_reads, fill = classification_type) + 
  geom_bar(stat = "identity", position = "stack") + 
  coord_flip() + 
  scale_fill_manual(values = cbPalette) + 
  theme_classic(base_size = 25) + 
  xlab("Sample") + 
  ylab("Number of Reads")

ggsave(snakemake@output[["added_reads_plot"]], height = 15, width = 20)

#### Stratified Heatmap of Relative Abundance Per Per Sample ####

conc_data_clean$species <- sub(" ", "_", conc_data_clean$species)

conc_data_clean_sep <- conc_data_clean %>%
  separate(species, into = c("Genus", "Species"), sep = "_") %>%
  filter(kraken_assigned_reads > filter)

conc_data_clean_genus_group <- conc_data_clean_sep %>%
  group_by(Genus) %>%
  summarise(n())

plot_stratified_heatmap <- function(){
  for(i in conc_data_clean_genus_group$Genus) {
  print(
    conc_data_clean_sep %>%
      filter(Genus == i) %>% 
      ggplot() + 
      aes(x = sample, y = Species) + 
      geom_tile(aes(fill = new_est_reads), linejoin = "round") + 
      scale_fill_gradientn(colours = c("white", "#F0E442", "#E69F00")) +
      facet_grid(rows = "Genus") + 
      theme_classic(base_size = 25) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Sample"))
  }
}

pdf(snakemake@output[["bracken_stratified_heatmap"]], height = 20, width = 20)
plot_stratified_heatmap() 
dev.off()

##### Heatmap of Species #####
plot_heatmap_per_sample <- function(data, filter_thresh){
  data %>%
    filter(new_est_reads > filter_thresh) %>%
    select(sample, species, fraction_total) %>%
    ggplot() +
    aes(x = sample, y = species) +
    geom_tile(aes(fill = log(fraction_total))) + 
    scale_fill_gradientn(colours = c("white", "orange", "blue")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + 
    ggtitle(paste("Log Relative Abundance of Microbial Species filter: Filtered > ", filter_thresh, sep = ""))
}

pdf(snakemake@output[["bracken_overall_heatmap"]], height = 20, width = 20)
plot_heatmap_per_sample(conc_data_clean, filter)
dev.off()
