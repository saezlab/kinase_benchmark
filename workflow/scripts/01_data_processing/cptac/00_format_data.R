#'

## Snakemake ---------------------------
if(exists("snakemake")){
  phosphosite_file <- snakemake@input$phospho
  pps_output <- snakemake@output$out
}else{
  phosphosite_file <- "data/datasets/cptac/brca_norm2prot_global_lm_log2_medCentRatio.rds"
  pps_output <- "results/01_processed_data/cptac/data/global/brca.csv"
}

## Libraries ---------------------------
library(tidyverse)

# Change naming and format to be consistent with normalised data
phosphosites <- readRDS(phosphosite_file)

df <- phosphosites %>%
  as.data.frame() %>%
  rownames_to_column("ID")

write_csv(df, pps_output)
