#!/usr/bin/env Rscript
library(tidyverse)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##### Import Data #####s

spp_tbl <- read_tsv(snakemake@input[["spp_table"]])
genus_tbl <- read_tsv(snakemake@input[["genus_table"]]) 
filter_threshold <- as.character(snakemake@params[["strat_thresh"]])
filter_threshold <- as.numeric(filter_threshold)

genus_tbl_clean <- genus_tbl %>%
  separate(study, into = c("sample", "junk"), sep = "_") %>%
  select(!junk)

spp_tbl_clean <- spp_tbl %>%
  separate(study, into = c("sample", "junk"), sep = "_") %>%
  select(!junk)
  
##### Plot Genus Composition #####

genus_thres <- as.character(snakemake@params[["genus_thresh"]])
genus_thres <- as.numeric(genus_thres)

genus_tbl_clean %>%
  filter(read_number > genus_thres) %>%
  ggplot() + 
  aes(x = sample, y = genus) + 
  geom_tile(aes(fill = log(read_number)), linejoin = "round") + 
  scale_fill_gradientn(colours = c("white", "orange", "blue")) +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(title = "Genus Composition (Threshold > 6000 reads)")

ggsave(snakemake@output[["genus_heatmap"]], height = 20, width = 15)

#### Plot Species Composition #### 
spp_thres <- as.character(snakemake@params[["species_thresh"]])
spp_thres <- as.numeric(spp_thres)

spp_tbl_clean %>%
  filter(read_number > spp_thres) %>%
  ggplot() + 
  aes(x = sample, y = species) + 
  geom_tile(aes(fill = log(read_number)), linejoin = "round") + 
  scale_fill_gradientn(colours = c("white", "orange", "blue")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(snakemake@output[["species_heatmap"]], height = 20, width = 15)

##### Faceted Plot By Genus #####
spp_tbl_split <- spp_tbl_clean %>% 
  separate(species, into = c("Genus", "Species1", "Species2"), sep = "_") %>%
  unite("Species", Species1:Species2, na.rm = TRUE) %>%
  filter(read_number > filter_threshold) 

genus_list <- spp_tbl_split %>%
  group_by(Genus) %>%
  summarise(n()) 

plot_all <- function() {
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
plot_all()
dev.off()

