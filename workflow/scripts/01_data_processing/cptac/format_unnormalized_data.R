#'

## Snakemake ---------------------------
if(exists("snakemake")){
  phosphosite_file <- snakemake@input$phospho
  pps_output <- snakemake@output$out
}else{
  phosphosite_file <- "data/CPTAC_original/brca_original_medcent_30plus.tsv"
  pps_output <- "data/CPTAC_phospho/final/brca_norm2prot_original_lm_log2_medCentRatio.rds"
}

## Libraries ---------------------------
library(tidyverse)

# Change naming and format to be consistent with normalised data
phosphosites <- read_tsv(phosphosite_file, col_types = cols())
saveRDS(phosphosites %>% column_to_rownames("site"), pps_output)
