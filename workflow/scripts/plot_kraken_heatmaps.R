#!/usr/bin/env Rscript
library(tidyverse)
library(MetBrewer)

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

##### Plot Function #####

kraken_heatmap <- function(clean_kraken_table, read_thresh, type){
  heatmap <- clean_kraken_table %>%
    filter(read_number > read_thresh) %>%
    filter(phylo_level != "Homo_sapiens") %>%
    ggplot() + 
    aes(x = sample, y = phylo_level) + 
    geom_tile(aes(fill = log(read_number)), linejoin = "round") + 
    scale_fill_gradientn(colours = met.brewer("Hiroshige", 4, type = "discrete")) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), 
          legend.position = "top", 
          legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15)) + 
    labs(fill = "Read Number (log)") + 
    xlab("Sample") + 
    ylab(type)
  
  return(heatmap)

}

species_heatmap <- kraken_heatmap(spp_tbl_clean, spp_thres, "Species")
genus_heatmap <- kraken_heatmap(genus_tbl_clean, genus_thres, "Genus")

ggsave(snakemake@output[["genus_heatmap"]], genus_heatmap, height = 20, width = 15)
ggsave(snakemake@output[["species_heatmap"]], species_heatmap, height = 20, width = 15)

##### Faceted Plot By Genus #####
spp_tbl_split <- spp_tbl_clean %>% 
  separate(phylo_level, into = c("Genus", "Species1", "Species2"), sep = "_") %>%
  unite("Species", Species1:Species2, na.rm = TRUE) 

genus_list <- spp_tbl_split %>%
  group_by(Genus) %>%
  summarise(n()) 

plot_species_facets <- function() {
  for(i in genus_list$Genus) {
    print(
      spp_tbl_split %>%
        filter(Genus == i) %>% 
        ggplot() + 
        aes(x = sample, y = Species) + 
        geom_tile(aes(fill = read_number), linejoin = "round") + 
        scale_fill_gradientn(colours = c("white", "#F0E442", "#E69F00")) +
        facet_grid(rows = "Genus") + 
        theme_classic(base_size = 25) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        )
    }
}

pdf(snakemake@output[["stratified_heatmap"]], height = 20, width = 20)
plot_species_facets()
dev.off()