### Load Packages
#!/usr/bin/env Rscript
library(tidyverse)

### Define and import data from {input} snakemake
columns <- c("Spp", "read_count")
data <- read_tsv(snakemake@input[["krak"]], col_names = columns)

### Hacky way to separate the Kraken out into R readable format
clean_output <- function(kraken_table){
  kraken_table_clean <- kraken_table %>%
    separate(Spp, 
             c("d", "d1", "domain", "p", "p1", "phylum", 
               "c", "c1", "class", "o", "o1", "order", "f", "f1", "familiy", "g", "g1", "genus")) %>%
    select(domain, phylum, class, order, familiy, genus, read_count) %>%
  
  return(kraken_table_clean)
}

### Do functions
data_clean <- clean_output(data) 

plot_krak <- function(data) {
  krak_plot <- function(data1) {
    ggplot(data1) + 
      aes(x = reorder(genus, read_count), y = read_count) +
      geom_bar(stat = "identity") + 
      coord_flip() + 
      theme_classic()
  } 
  
  sum_data <- data %>%
    summarise(mean_read_count = mean(read_count))
  
  if(sum_data$mean_read_count > 3000) {
    data %>%
      filter(read_count > 50000) %>%
      select(genus, read_count) %>%
      krak_plot()
  } else if(sum_data$mean_read_count < 3000) {
    data %>%
      filter(read_count > 4000) %>%
      select(genus, read_count) %>%
      krak_plot()
  }

}

pdf(snakemake@output[["plot"]])
plot_krak(data_clean)
dev.off()
