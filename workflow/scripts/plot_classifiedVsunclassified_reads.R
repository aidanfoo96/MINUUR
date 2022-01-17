#!/usr/bin/env Rscript
library(tidyverse)
library(MetBrewer)

#### Import Kraken Data
concatenated_kraken_summaries <- read_tsv(snakemake@input[["concatenated_kraken_summary"]], 
                                            col_names = c("proportion_of_reads", 
                                                        "number_of_reads", 
                                                        "number_of_reads2",
                                                        "read_code", 
                                                        "read_code2", 
                                                        "classification_status", 
                                                        "file_path"))

concatenated_kraken_summaries_clean <- concatenated_kraken_summaries %>%
  separate(file_path, into = c("junk", "filename"), sep = "classified_summary/") %>%
  separate(filename, into = c("sample", "junk2"), sep = "_classified_summary.txt") %>%
  select(sample, proportion_of_reads, number_of_reads, classification_status) 


#### Plot Proportion of Classified vs Unclassified Reads from Kraken
plot_kraken_proportion <- function(kraken_table){
  
  concatenated_kraken_summaries_clean_plot <- kraken_table %>%
    ggplot() + 
    aes(x = sample, y = number_of_reads, fill = classification_status) + 
    geom_bar(stat = "identity", position = "fill") + 
    theme_bw(base_size = 15) + 
    scale_fill_manual(values = met.brewer("Redon", n = 2, type = "continuous"), 
                      labels = c("Classified Reads", "Unclassified Reads")) + 
    xlab("Sample") + 
    ylab("Proportion of Reads") + 
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), 
          legend.position = "top")
  
  return(concatenated_kraken_summaries_clean_plot)
  
}

kraken_classification_proportions <- plot_kraken_proportion(concatenated_kraken_summaries_clean)

#### Save file
ggsave(snakemake@output[["classified_proportions"]], kraken_classification_proportions)
