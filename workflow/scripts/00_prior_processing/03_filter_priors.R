#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$raw
  kinhub_file <- snakemake@input$kin
  output_file <- snakemake@output$filter
}else{
  input_file <- "results/00_prior/raw/phosformer5.tsv"
  kinhub_file <- "data/misc/kinase_list.tsv"
  output_file <- "results/00_prior/phosformer5.tsv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

## Load data ---------------------------
net <- read_tsv(input_file, col_types = cols())
kin_list <- read_tsv(kinhub_file, col_types = cols())

kinases <- kin_list$kinase

## Check kinases not covered in data ---------------------------
net %>%
  filter(!source %in% kinases) %>%
  pull(source) %>%
  unique()

net_filtered <- net %>% filter(source %in% kinases)

write_tsv(net_filtered, output_file)
