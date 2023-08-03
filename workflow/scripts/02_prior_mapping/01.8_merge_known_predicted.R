if(exists("snakemake")){
  known_file <- snakemake@input$known_file
  predicted_file <- snakemake@input$predicted_file
  known_decryptm <- snakemake@input$known_decryptm
  predicted_dectyptm <- snakemake@input$predicted_dectyptm
  output_file <- snakemake@output$tsv
  output_decryptm <- snakemake@output$tsv_decryptm
}else{
  known_file <- "results/prior/GPS.tsv"
  predicted_file <- "results/prior/networkin.tsv"
  known_decryptm <- "results/decryptm/prior/GPS.tsv"
  predicted_dectyptm <- "results/decryptm/prior/networkin.tsv"
  output_file <- "results/prior/GPS_networkin.tsv"
  output_decryptm <- "results/decryptm/prior/GSP_networkin.tsv"
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
merge_df <- merge_df[!duplicated(merge_df),]

write_tsv(merge_df, output_file)

# Decryptm
known_dc <- read_tsv(known_decryptm)
predicted_dc <- read_tsv(predicted_dectyptm)

## only add pps from kinases already present in known targets
predicted_dc_kin <- predicted_dc %>%
  dplyr::filter(source %in% known_dc$source)

merge_df_decryptm <- rbind(known_dc, predicted_dc_kin)
merge_df_decryptm <- merge_df_decryptm[!duplicated(merge_df_decryptm),]

write_tsv(merge_df_decryptm, output_decryptm)
