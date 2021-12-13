#!/usr/bin/env Rscript
library(tidyverse)

metaphlan_tbl <- read_tsv(snakemake@input[["concat_tbl"]], skip = 1) %>%
  select(!NCBI_tax_id)

#### Get Kingdom ####
metaphlan_tbl_king <- metaphlan_tbl %>%
  filter(!grepl("p__|g__|s__|c__|o__|f__", clade_name)) %>%
  mutate_at("clade_name", str_replace, "k__", "")

metaphlan_tbl_king_long <- metaphlan_tbl_king %>%
  pivot_longer(cols = !clade_name, 
               names_to = "samples", 
               values_to = "rel_abund")

#### Get Genus ####
metaphlan_tbl_genus <- metaphlan_tbl %>%
  filter(grepl("g__", clade_name)) %>%
  separate(clade_name, c("junk", "genus"), sep = "g__") %>%
  select(!junk) %>%
  filter(!grepl("s__", genus))

metaphlan_tbl_genus_long <- metaphlan_tbl_genus %>% 
  pivot_longer(cols = !genus, 
               names_to = "samples", 
               values_to = "rel_abund")

#### Get Species ####
metaphlan_tbl_spp <- metaphlan_tbl %>%
  filter(grepl("s__", clade_name)) %>%
  separate(clade_name, c("junk", "species"), sep = "s__") %>%
  select(!junk)

metaphlan_tbl_spp_long <- metaphlan_tbl_spp %>% 
  pivot_longer(cols = !species,
               names_to = "samples", 
               values_to = "rel_abund")

#### Export Dataframes to TSV files ####
write.table(metaphlan_tbl_king_long, file = snakemake@output[["kingdom_table"]], sep = "\t", row.names = FALSE)
write.table(metaphlan_tbl_genus_long, file = snakemake@output[["genus_table"]], sep = "\t", row.names = FALSE)
write.table(metaphlan_tbl_spp_long, file = snakemake@output[["spp_table"]], sep = "\t", row.names = FALSE)
