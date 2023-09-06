if(exists("snakemake")){
  gps_file <- snakemake@input$gps
  ppsp_file <- snakemake@input$ppsp
  gps_decryptm <- snakemake@input$gps_decryptm
  ppsps_decryptm <- snakemake@input$ppsp_decryptm
  output_file <- snakemake@output$tsv
  output_decryptm <- snakemake@output$tsv_decryptm
}else{
  gps_file <- "results/prior/GPS.tsv"
  ppsp_file <- "results/prior/phosphositeplus.tsv"
  gps_decryptm <- "results/decryptm/prior/GPS.tsv"
  ppsp_dectyptm <- "results/decryptm/prior/phosphositeplus.tsv"
  output_file <- "results/prior/GPS_PPSP.tsv"
  output_decryptm <- "results/decryptm/prior/GSP_PPSP.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## load and merge files ---------------------------
gps <- read_tsv(gps_file)
ppsp <- read_tsv(ppsp_file)

merge_df <- rbind(gps, ppsp)
merge_df <- merge_df[!duplicated(merge_df),]

write_tsv(merge_df, output_file)

gps_dc <- read_tsv(gps_decryptm)
ppsp_dc <- read_tsv(ppsp_file)

merge_df_decryptm <- rbind(gps_dc, ppsp_dc)
merge_df_decryptm <- merge_df_decryptm[!duplicated(merge_df_decryptm),]

write_tsv(merge_df_decryptm, output_decryptm)
