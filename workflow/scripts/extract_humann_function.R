### Load Packages
#!/usr/bin/env Rscript
library(tidyverse)

genes <- read_tsv(snakemake@input[["renamed_table"]])
role <- as.character(snakemake@params[["role"]])

genes_clean <- genes %>%
  separate(`# Gene Family`, c("gene", "genus"), sep = "g__") %>%
  separate(genus, c("genus", "species"), sep = "s__") %>%
  separate(gene, c("label", "description"), sep = ":") %>%
  select(!genus) 

genes_clean_long <- genes_clean %>%
  pivot_longer(cols = !c(label, description, species), 
               names_to = "sample", 
               values_to = "rel_abund")


genes_clean_long %>% 
  filter(grepl(role, description)) %>%
  filter(rel_abund > 0) %>%
  ggplot() + 
  aes(x = description, fill = species) %>%
  geom_bar(position = "stack") + 
  coord_flip() + 
  theme_minimal()

ggsave(snakemake@output[["stacked_bar_plot"]], height = 20, width = 20)

genes_clean_long %>%
  filter(rel_abund > 0) %>%
  group_by(species, description) %>%
  summarise(sum = n()) %>%
  filter(grepl(role, description)) %>%
  ggplot() + 
  aes(x = species, y = description) + 
  geom_tile(aes(fill = sum), colour = "black") + 
  scale_fill_gradientn(colours = c("yellow", "orange", "blue")) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(snakemake@output[["heatmap_plot"]], height = 20, width = 20)
