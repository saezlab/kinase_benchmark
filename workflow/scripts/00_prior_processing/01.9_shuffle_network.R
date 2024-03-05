#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$tsv
  output_file <- snakemake@output$shuffled
}else{
  input_file <- "results/prior/raw/phosphositeplus.tsv"
  output_file <- "results/prior/raw/shuffled.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## Shuffle network ---------------------------
set.seed(1234)

net <- read_tsv(input_file,
                col_types = cols())

net_shuffled <- transform(net, source = sample(source))

## Save shuffled network ---------------------------
write_tsv(net_shuffled, output_file)
