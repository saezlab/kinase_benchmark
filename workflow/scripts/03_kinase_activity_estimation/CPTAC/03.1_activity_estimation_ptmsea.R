if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  dataset_name <- snakemake@wildcards$dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  normalisation <- snakemake@wildcards$normalisation
  log <- snakemake@output$rds
  output_folder <- snakemake@params$output_folder
}else{
  dataset <- "results/datasets/gbm.gct"
  dataset_name <- "gbm"
  PKN <- "results/prior/ptm-sea/GPS.gmt"
  PKN_name <- "GPS"
  output_folder <- "results/activity_scores_ptmsea"
  log <- "results/activity_scores_ptmsea/log/gbm_GPS.log"
  if(!require("ssGSEA2")) remotes::install_github('nicolerg/ssGSEA2', upgrade='never')
}

## Libraries ---------------------------
library(ssGSEA2)
library(tidyselect)

## Run ssGSEA/PTM-SEA ---------------------------
res <- run_ssGSEA2(dataset,
                  output.prefix = paste(normalisation, dataset_name, PKN_name, sep = "_"),
                  gene.set.databases = PKN,
                  output.directory = output_folder,
                  sample.norm.type = "none",
                  weight = 0.75,
                  correl.type = "rank",
                  statistic = "area.under.RES",
                  output.score.type = "NES",
                  nperm = 50,
                  min.overlap = 5,
                  extended.output = TRUE,
                  global.fdr = FALSE,
                  log.file = log)
