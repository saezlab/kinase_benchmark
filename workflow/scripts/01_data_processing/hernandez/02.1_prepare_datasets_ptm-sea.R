if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$gct
}else{
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
  output_file <- "results/hernandez/datasets/benchmark.gct"
}

## Libraries ---------------------------
library(cmapR)
library(tidyverse)

## Prepare files for PTM-SEA ---------------------------
dataset <- read_csv(file = dataset, col_types = cols())

dataset <- dataset %>%
  column_to_rownames("ID")
dataset_m <- as.matrix(dataset)

dataset_gct <- new("GCT", mat=dataset_m)

## Save gct ---------------------------
write_gct(dataset_gct, output_file, appenddim = FALSE)
