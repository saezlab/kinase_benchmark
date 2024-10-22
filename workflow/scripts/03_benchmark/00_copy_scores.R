
## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  input_hernandez <- snakemake@input$hernandez
  output_file <- snakemake@output$output
}else{
  input_file <- "results/02_activity_scores/hijazi/scores/GPS.rds"
  output_file <- "results/02_activity_scores/hijaziDiscoverX/scores/GPS.rds"
}

## Libraries ---------------------------
library(tidyverse)

## Load scores and meta ---------------------------
act_scores_hijazi <- readRDS(input_file)

## Save merged data ---------------------------
saveRDS(act_scores_hijazi, output_file)
