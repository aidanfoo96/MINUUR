#!/usr/bin/env Rscript
library(tidyverse)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

data <- read_tsv(snakemake@input[["long_read_tbl"]])

plot <- data %>% 
  filter(sequence_type == c("reads unmapped:", "classified_reads")) %>%
  ggplot() + 
  aes(x = sample, y = num_reads, fill = sequence_type) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_brewer(palette="Set2", labels = c("Classified Reads", "Unclassified Reads")) + 
  coord_flip() + 
  ylab("Proportion of Reads") + xlab("Sample") + 
  labs(title = "Proportion of Classified to Unclassified Reads") + 
  guides(fill = guide_legend("Sequence Type")) + 
  theme_minimal()
  
ggsave(snakemake@output[["classified_proportions"]], height = 15, width = 15)  