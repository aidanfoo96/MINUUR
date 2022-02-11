### Load Packages
#!/usr/bin/env Rscript
library(tidyverse)
library(MetBrewer)

#### Load Data ####
concatenated_alignments <- read_tsv(snakemake@input[["alignment_stats"]], col_names = c("SN", "sequence", "num_sequences", "file_samples"))
concatenated_alignments$file_samples[concatenated_alignments$file_samples == "# excluding supplementary and secondary reads"] <- ""

#### Get List of Sample Names ####
sample_names <- concatenated_alignments %>%
  separate(file_samples, c("file", "sample"), sep = "stats_ffq/") %>%
  separate(sample, c("sample", "junk"), sep = "_align") %>%
  select(sequence, num_sequences, sample) %>%
  group_by(sample) %>%
  summarise(n()) %>%
  select(sample) %>%
  filter(sample != "NA")

sample_names <- as.vector(sample_names$sample)

#### Clean Data ####
concatenated_alignments_clean <- concatenated_alignments %>%
  select(!SN) %>%
  separate(file_samples, c("file", "sample"), sep = "stats_ffq/") %>%
  separate(sample, c("sample", "junk"), sep = "_align") %>%
  select(sequence, num_sequences, sample) %>%
  mutate(sample = replace(sample, is.na(sample), sample_names))

##### Generate Summary Tables ######
# Wider, 'human readable' format with proportions
concatenated_alignments_clean_wide <- concatenated_alignments_clean %>% 
  group_by(sample) %>%
  pivot_wider(names_from = sequence, values_from = num_sequences) %>%
  mutate(proportion_mapped = `reads mapped:` / `raw total sequences:` * 100) 

write.table(concatenated_alignments_clean_wide, file = snakemake@output[["human_read_tbl"]], sep = "\t", row.names = FALSE)
write.table(concatenated_alignments_clean, file = snakemake@output[["long_read_tbl"]], sep = "\t", row.names = FALSE)

#### Plot Summaries ####
alignment_statistics_plot <- concatenated_alignments_clean %>%
  ggplot() + 
  aes(x = sample, y = num_sequences, fill = sequence) %>%
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw(base_size = 15) + 
  scale_fill_manual(values = met.brewer("Redon", n = 3, type = "discrete"), 
                    labels = c("Raw Read Number", "Reads Mapped", "Reads Unmapped")) + 
  xlab("Sample") + 
  ylab("Read Number") + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), 
        legend.position = "top") 

ggsave(snakemake@output[["alignment_stat_plot"]], alignment_statistics_plot)