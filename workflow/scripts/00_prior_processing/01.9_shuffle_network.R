#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$tsv
  output_file <- snakemake@output$shuffled
}else{
  input_file <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  output_file <- "results/01_processed_data/hernandez/mapped_priors/shuffled.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## Shuffle network ---------------------------
set.seed(123)

net <- read_tsv(input_file,
                col_types = cols())

net_shuffled <- transform(net, source = sample(source)) %>%
  distinct()

## Save shuffled network ---------------------------
write_tsv(net_shuffled, output_file)
