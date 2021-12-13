### Load Packages
#!/usr/bin/env Rscript
library(tidyverse)

#### Load Data ####
krak_class <- read_tsv(snakemake@input[["classified_reads"]])
alignment_sum <- read_tsv(snakemake@input[["alignment_stats"]], col_names = c("SN", "sequence_type", "no.reads", "sample"))

#### Clean Data ####
alignment_sum_clean <- alignment_sum %>%
  separate(sample, c("junkn", "junk2", "junk3", "sample"), sep = "/") %>%
  separate(sample, c("sample", "junk"), sep = "_") %>%
  select(!c(junkn, junk2, junk3, junk, SN)) %>%
  pivot_wider(names_from = "sequence_type", 
              values_from = "no.reads")
  

krak_class_clean <- krak_class %>%
  separate(study, c("sample", "junk"), sep = "_") %>%
  select(!(junk)) %>%
  rename(classified_reads = "total_reads")
  
#### Join Kraken Results and Alignment Summaries ####
krak_alignments_join <- alignment_sum_clean %>%
  left_join(krak_class_clean)

##### Generate Summary Tables ######
# Table suitable for ggplot - long format
krak_alignments_join_long <- krak_alignments_join %>% 
  pivot_longer(cols = c(`raw total sequences:`, `reads mapped:`,
                        `reads unmapped:`, classified_reads), 
               names_to = "sequence_type", 
               values_to = "num_reads")

# Wider, 'human readable' format with proportions
krak_alignments_join_clean <- krak_alignments_join %>%
  rename(raw_sequence_num = `raw total sequences:`, 
         mapped_reads = `reads mapped:`,
         unmapped_reads = `reads unmapped:`) %>%
  mutate(unclassified_reads = (unmapped_reads - classified_reads)) %>%
  mutate(proportion_unmapped = (classified_reads/unmapped_reads)*100)

write.table(krak_alignments_join_clean, file = snakemake@output[["human_read_tbl"]], sep = "\t", row.names = FALSE)
write.table(krak_alignments_join_long, file = snakemake@output[["long_read_tbl"]], sep = "\t", row.names = FALSE)
