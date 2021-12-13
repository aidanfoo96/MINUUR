#!/usr/bin/env Rscript
library(tidyverse)

genus_delim <- function(){
  if(snakemake@params[["kraken_db_type"]]){
    return("g__")
  }else{
    return("g__g__")
  }
}
  

spp_delim <- function(){
  if(snakemake@params[["kraken_db_type"]]){
    return("s__")
  }else{
    return("s__s__")
  }
}
  
genus <- genus_delim()
spp <- spp_delim()

data <- read_tsv(snakemake@input[["combined_report"]])

##### Get Domain #####
kingdom_data <- data %>%
  filter(!grepl("p__|g__|s__|c__|o__|f__", `#Classification`)) %>%
  separate(`#Classification`, c("junk", "kingdom"), sep = "^k__") %>%
  select(!junk)

kingdom_data_long <- kingdom_data %>%
  pivot_longer(cols = !kingdom, names_to = "study", values_to = "read_number") 

##### Get genus  #####
genus_data <- data %>%
  filter(grepl("g__", `#Classification`)) %>%
  separate(`#Classification`, c("junk", "genus"), sep = genus) %>%
  select(!junk) %>%
  filter(!grepl("s__", genus)) 

genus_data_long <- genus_data %>%
  pivot_longer(cols = !genus, names_to = "study", values_to = "read_number") 

##### Get species  #####
species_data <- data %>% 
  filter(grepl("s__", `#Classification`)) %>%
  separate(`#Classification`, c("junk", "species"), sep = spp) %>%
  select(!junk)

species_data_long <- species_data %>% 
  pivot_longer(cols = !species, names_to = "study", values_to = "read_number") 

##### How many reads per study? #####
total_reads_per_study <- species_data_long %>% 
  group_by(study) %>%
  summarise(total_reads = sum(read_number)) 

classified_reads_plot <- total_reads_per_study %>% 
  ggplot() + 
  aes(x = study, y = total_reads) + 
  geom_col() + 
  theme_classic() + 
  coord_flip()

##### Species Classification Passed a specific threshold (unfinished) ##### 
species_data_filtered <- species_data_long %>%
  filter(read_number > 1000)

genus_data_filtered <- genus_data_long %>%
  filter(read_number > 1000)

#### Export Dataframes to TSV files ####
write.table(kingdom_data_long, file = snakemake@output[["kingdom_table"]], sep = "\t", row.names = FALSE)
write.table(genus_data_long, file = snakemake@output[["genus_table"]], sep = "\t", row.names = FALSE)
write.table(species_data_long, file = snakemake@output[["spp_table"]], sep = "\t", row.names = FALSE)
write.table(total_reads_per_study, file = snakemake@output[["classified_reads"]], sep = "\t", row.names = FALSE)

#### Export Classified Reads Plot ####
pdf(snakemake@output[["classified_reads_plot"]])
print(classified_reads_plot)
dev.off()

