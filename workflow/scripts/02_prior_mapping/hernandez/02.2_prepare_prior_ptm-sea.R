if(exists("snakemake")){
  PKN <- snakemake@input$file_PKN
  output_file <- snakemake@output$gmt
}else{
  PKN <- "results/hernandez/prior/iKiPdb.tsv"
  output_file <- "results/hernandez/prior/ptm-sea/iKiPdb.gmt"
}


## Libraries ---------------------------
library(rWikiPathways)
library(tidyverse)


## Prepare prior files for PTM-SEA ---------------------------
prior <- read_tsv(PKN, col_types = cols()) %>%
  dplyr::select(source, target) %>%
  as.data.frame()

## Save gmt ---------------------------
rWikiPathways::writeGMT(prior, output_file)
