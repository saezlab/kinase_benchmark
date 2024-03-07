#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$scores
  subset <- snakemake@input$filter
  meta_file <- snakemake@input$meta
  output_file <- snakemake@output$output
  output_meta <- snakemake@output$meta_out
}else{
  input_file <- "results/03_benchmark/merged/01_input_bench/fgsea-GPS.csv"
  subset <- "results/03_benchmark/merged/01_input_bench_subset//tmp/filter_subset.csv"
  meta_file <- "results/03_benchmark/merged/01_input_bench/obs_fgsea-GPS.csv"
  output_file <- "results/03_benchmark/merged/01_input_bench/tmp/fgsea-GPS.csv"
  output_meta <- "results/03_benchmark/merged/01_input_bench/tmp/obs_fgsea-GPS.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Load data ---------------------------
df <- read_csv(input_file, col_types = cols())
meta <- read_csv(meta_file, col_types = cols())
shared <- read_csv(subset, col_types = cols()) %>%
  pull(id)

## Generate shared subset ---------------------------
subset_df <- df %>%
  pivot_longer(!experiment, values_to = "score", names_to = "kinase") %>%
  mutate(id = paste(experiment, kinase, sep = ":")) %>%
  filter(id %in% shared) %>%
  dplyr::select(-id) %>%
  pivot_wider(values_from = "score", names_from = "kinase")

write_csv(subset_df, output_file)
write_csv(meta, output_meta)
