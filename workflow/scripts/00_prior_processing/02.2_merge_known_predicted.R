if(exists("snakemake")){
  known_file <- snakemake@input$known_file
  predicted_file <- snakemake@input$predicted_file
  output_file <- snakemake@output$tsv
}else{
  known_file <- "results/00_prior/phosphositeplus.tsv"
  predicted_file <- "results/00_prior/iKiPdb.tsv"
  output_file <- "results/00_prior/merged/phosphositeplus_iKiPdb.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## load table for kinase targets (GPS Gold Standard (GS) set) ---------------------------
known <- read_tsv(known_file, col_types = cols())
predicted <- read_tsv(predicted_file, col_types = cols())

## Construct kinase-substrate interaction network ---------------------------
predicted <- predicted %>%
  filter(source %in% known$source) # only add targets for existing kinases

merged <- rbind(known, predicted) %>%
  distinct(source, target, .keep_all = T)

write_tsv(merged, output_file)
