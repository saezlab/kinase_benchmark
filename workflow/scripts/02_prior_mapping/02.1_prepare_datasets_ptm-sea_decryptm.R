if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$gct
}else{
  dataset <- "results/decryptm/processed_data/R2_pEC50.csv"
  output_file <- "results/decryptm/datasets/decryptm.gct"
}

## Libraries ---------------------------
library(cmapR)
library(tidyverse)

## Prepare files for PTM-SEA ---------------------------
dataset <- read_tsv(file = dataset, col_types = cols())

dataset <- dataset %>%
  column_to_rownames("pps_id")
dataset_m <- as.matrix(dataset)

dataset_gct <- new("GCT", mat=dataset_m)

## Save gct ---------------------------
write_gct(dataset_gct, output_file, appenddim = FALSE)
