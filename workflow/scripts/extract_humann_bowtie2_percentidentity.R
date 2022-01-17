### Load Packages
#!/usr/bin/env Rscript
library(tidyverse)

krak_class <- read_tsv(snakemake@input[["classified_reads"]])

bowtie2_aligned <- read_tsv(snakemake@input[["bowtie2_humann_align"]], 
                            col_names = c("query_accession", "reference_sequence_name", 
                                          "percentage_identity", "alignment_length", 
                                          "empty1", "empty2", 
                                          "empty3", "empty4", 
                                          "empty5", "empty6", 
                                          "e_value","empty7"))

sample_name <- snakemake@input[["bowtie2_humann_align"]]

bowtie2_aligned_sel <- bowtie2_aligned %>% 
  select(c(query_accession, 
           reference_sequence_name, 
           percentage_identity, 
           alignment_length,
           e_value))

percent_iden_grouped <- bowtie2_aligned_sel %>%
  group_by(percentage_identity) %>%
  summarise(gene_number = n()) %>%
  mutate(sample_name = sample_name) %>%
  separate(sample_name, into = c("sample", "junk"), sep = "_unmapped") %>%
  select(!junk)

write.table(percent_iden_grouped, file = snakemake@output[["summarised_genes"]], sep = "\t", row.names = FALSE)
