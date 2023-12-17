if(exists("snakemake")){
  known_file <- snakemake@input$known_file
  predicted_file <- snakemake@input$predicted_file
  output_file <- snakemake@output$tsv
}else{
  known_file <- "results/prior/raw/GPS.tsv"
  predicted_file <- "results/prior/raw/networkin.tsv"
  output_file <- "results/prior/raw/GPS_networkin.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## load table for kinase targets (GPS Gold Standard (GS) set) ---------------------------
known <- read_tsv(known_file, col_types = cols())
predicted <- read_tsv(predicted_file, col_types = cols())

## Construct kinase-substrate interaction network ---------------------------
merged <- rbind(known, predicted) %>%
  distinct(source, target, .keep_all = T)

write_tsv(merged, output_file)
