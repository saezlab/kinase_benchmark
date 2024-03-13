#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input
  output_file <- snakemake@output
}else{
  input_file <- "results/00_prior/phosphositeplus.tsv"
  net_file <- "results/00_prior/merged/phosphositeplus_networkin.tsv"
  ikip_file <- "results/00_prior/merged/phosphositeplus_iKiPdb.tsv"
  meta_file <- "results/01_processed_data/merged/data/benchmark_metadata.csv"
  output_file <- "path_to_output"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
meta <- read_csv(meta_file, col_types = cols())

phosphositeplus <- read_tsv(input_file, col_types = cols())
net_merged <- read_tsv(net_file, col_types = cols())
ikip_merged <- read_tsv(ikip_file, col_types = cols())

n_ppsp <- phosphositeplus %>%
  group_by(source) %>%
  summarise(phosphositeplus = n())

n_net <- net_merged %>%
  group_by(source) %>%
  summarise(networkin = n())

n_ikip <- ikip_merged %>%
  group_by(source) %>%
  summarise(ikip = n())

added_targets <- reduce(list(n_ppsp, n_net, n_ikip), full_join) %>%
  filter(!is.na(phosphositeplus)) %>%
  filter(phosphositeplus >= 3) %>%
  mutate(networkin = networkin/phosphositeplus) %>%
  mutate(ikipdb = ikip/phosphositeplus) %>%
  dplyr::select(source, networkin, ikipdb) %>%
  pivot_longer(!source, values_to = "percentage", names_to = "net")

added_targets <- added_targets %>%
  mutate(bench = case_when(
    source %in% meta$target ~ "bench",
    !source %in% meta$target ~ "background"
  ))

ggplot(added_targets, aes(x = net, y = percentage, fill = bench)) +
  geom_boxplot() +
  ylab("ratio")

