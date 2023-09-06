if(exists("snakemake")){
  known_file <- snakemake@input$known_file
  predicted_file <- snakemake@input$predicted_file
  output_file <- snakemake@output$tsv
}else{
  known_file <- "results/hernandez/prior/GPSppsp.tsv"
  predicted_file <- "results/hernandez/prior/networkin.tsv"
  output_file <- "results/hernandez/prior/GPSppsp_networkin.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## load and merge files ---------------------------
known <- read_tsv(known_file)
predicted <- read_tsv(predicted_file)

## only add pps from kinases already present in known targets
predicted_kin <- predicted %>%
  dplyr::filter(source %in% known$source)

merge_df <- rbind(known, predicted_kin)
merge_df <- merge_df %>%
  distinct()

write_tsv(merge_df, output_file)
