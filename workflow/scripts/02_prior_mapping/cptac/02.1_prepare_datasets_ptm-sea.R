if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$gct
}else{
  dataset <- "data/CPTAC_phospho/hnscc_phospho_data_median_centered.tsv"
  output_file <- "results/datasets/hnscc.gct"
}

## Libraries ---------------------------
library(cmapR)
library(tidyverse)

## Prepare files for PTM-SEA ---------------------------
dataset <- readRDS(file)

dataset_m <- as.matrix(dataset)

dataset_gct <- new("GCT", mat=dataset_m)

## Save gct ---------------------------
write_gct(dataset_gct, output_file, appenddim = FALSE)
