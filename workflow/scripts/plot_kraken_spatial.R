#!/usr/bin/env Rscript
library(tidyverse)
library(MetBrewer)
library(treemapify)

##### Import Data #####

spp_tbl <- read_tsv(snakemake@input[["spp_table"]])
genus_tbl <- read_tsv(snakemake@input[["genus_table"]]) 

##### Define User Read Threshold #####
filter_threshold <- as.character(snakemake@params[["strat_thresh"]])
filter_threshold <- as.numeric(filter_threshold)

genus_thres <- as.character(snakemake@params[["genus_thresh"]])
genus_thres <- as.numeric(genus_thres)

spp_thres <- as.character(snakemake@params[["species_thresh"]])
spp_thres <- as.numeric(spp_thres)

##### Clean and Prep Data #####
genus_tbl_clean <- genus_tbl %>%
  separate(study, into = c("sample", "junk"), sep = "_") %>%
  select(!junk) %>%
  rename("phylo_level" = genus)

spp_tbl_clean <- spp_tbl %>%
  separate(study, into = c("sample", "junk"), sep = "_") %>%
  select(!junk) %>%
  rename("phylo_level" = species)

##### Plot Spatial Plots #####
plot_spatial_classification <- function(classificaiton_tbl, read_thresh){
  tbl_clean_sum <- classificaiton_tbl %>%
    group_by(phylo_level) %>%
    summarise(summed_read_num = sum(read_number)) %>%
    filter(summed_read_num > read_thresh)
  
  spatial_plot <- tbl_clean_sum %>%
    filter(phylo_level != "Homo_sapiens") %>%
    filter(summed_read_num > read_thresh) %>%
    ggplot(aes(area = summed_read_num, fill = summed_read_num, label = phylo_level, subgroup = phylo_level)) +
    geom_treemap() + 
    geom_treemap_subgroup_border(colour = "white", size = 1) + 
    geom_treemap_text(colour = "white", place = "center", reflow = T) + 
    scale_fill_gradientn(colours = met.brewer("Hiroshige", 10, type = "discrete"), name = "Total Read Number") +
    theme_minimal(base_size = 15) + 
    theme(legend.position = "right")
  
  return(spatial_plot)
}

species_spatial_plot <- plot_spatial_classification(spp_tbl_clean, spp_thres)
genus_spatial_plot <- plot_spatial_classification(genus_tbl_clean, genus_thres)

ggsave(snakemake@output[["genus_spatial_plot"]], genus_spatial_plot, height = 20, width = 15)
ggsave(snakemake@output[["species_spatial_plot"]], species_spatial_plot, height = 20, width = 15)