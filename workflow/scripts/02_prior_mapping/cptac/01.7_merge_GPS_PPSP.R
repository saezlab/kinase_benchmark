if(exists("snakemake")){
  gps_file <- snakemake@input$gps
  ppsp_file <- snakemake@input$ppsp
  output_file <- snakemake@output$tsv
}else{
  gps_file <- "results/cptac/prior/GPS.tsv"
  ppsp_file <- "results/cptac/prior/phosphositeplus.tsv"
  output_file <- "results/cptac/prior/GPS_PPSP.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## load and merge files ---------------------------
gps <- read_tsv(gps_file)
ppsp <- read_tsv(ppsp_file)

merge_df <- rbind(gps, ppsp)
merge_df <- merge_df[!duplicated(merge_df),]

write_tsv(merge_df, output_file)
