
## Snakemake ---------------------------
if(exists("snakemake")){
  meta_file <- snakemake@input$meta
  meta_hernandez <- snakemake@input$metah
  meta_out <- snakemake@output$meta_out
}else{
  meta_file <- "results/01_processed_data/hijazi/data/benchmark_metadata.csv"
  meta_hernandez <- "results/01_processed_data/hernandez/data/benchmark_metadata.csv"
  meta_out <- "results/01_processed_data/merged/data/benchmark_metadata.csv"
}

## Libraries ---------------------------
library(tidyverse)

obs_hijazi <- read_csv(meta_file, col_types = cols()) %>%
  dplyr::select(id, sign, target)
obs_hernandez <- read_csv(meta_hernandez, col_types = cols()) %>%
  dplyr::select(id, sign, target)
obs <- rbind(obs_hijazi, obs_hernandez)

obs <- separate_rows(obs, target, sep = ";")

write_csv(obs, meta_out)
