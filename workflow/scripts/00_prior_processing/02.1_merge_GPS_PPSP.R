if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  GPS_file <- snakemake@input$gps
  ptmsig_file <- snakemake@input$ptmsig
  output_file_merged <- snakemake@output$merged
}else{
  ppsp_file <- "results/prior/raw/phosphositeplus.tsv"
  GPS_file <- "results/prior/raw/GPS.tsv"
  ptmsig_file <- "results/prior/raw/ptmsigdb.tsv"
  output_file_merged <- "results/prior/raw/GSknown.tsv"
}

## Libraries ---------------------------
library(tidyverse)

## load table for kinase targets (GPS Gold Standard (GS) set) ---------------------------
pps <- read_tsv(ppsp_file, col_types = cols())
gps <- read_tsv(GPS_file, col_types = cols())
ptm <- read_tsv(ptmsig_file, col_types = cols())

## Construct kinase-substrate interaction network ---------------------------
GPS_ppsp <- rbind(pps, gps, ptm) %>%
  distinct(source, target, .keep_all = T)

write_tsv(GPS_ppsp, output_file_merged)
