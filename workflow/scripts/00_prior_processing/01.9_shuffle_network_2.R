#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$tsv
  output_file <- snakemake@output$shuffled
}else{
  input_file <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  output_file <- "results/01_processed_data/hernandez/mapped_priors/shuffled2.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## Shuffle network ---------------------------
set.seed(123)

net <- read_tsv(input_file,
                col_types = cols())

net_m <- net %>%
  pivot_wider(names_from = target, values_from = mor) %>%
  column_to_rownames("source")

colnames(net_m) <- sample(colnames(net_m))

net_shuffled <- net_m %>%
  rownames_to_column("source") %>%
  pivot_longer(-source, names_to = "target", values_to = "mor") %>%
  filter(!is.na(mor))

## Save shuffled network ---------------------------
write_tsv(net_shuffled, output_file)
