if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  log <- snakemake@output$rds
  output_folder <- snakemake@params$output_folder
  minsize <- snakemake@params$minsize
}else{
  dataset <- "results/01_processed_data/hernandez/datasets/benchmark.gct"
  PKN <- "results/01_processed_data/hernandez/mapped_priors/ptm-sea/phosphositeplus.gmt"
  PKN_name <- "phosphositeplus"
  output_folder <- "results/02_activity_scores/hernandez/ptmsea"
  log <- "results/02_activity_scores/hernandez/ptmsea/log/phosphositeplus.log"
  if(!require("ssGSEA2")) remotes::install_github('nicolerg/ssGSEA2', upgrade='never')
  minsize <- 1
}
minsize <- as.double(minsize)

## Libraries ---------------------------
library(ssGSEA2)
library(tidyselect)

## Run ssGSEA/PTM-SEA ---------------------------
res <- run_ssGSEA2(dataset,
                  output.prefix = PKN_name,
                  gene.set.databases = PKN,
                  output.directory = output_folder,
                  sample.norm.type = "none",
                  weight = 0.75,
                  correl.type = "rank",
                  statistic = "area.under.RES",
                  output.score.type = "NES",
                  nperm = 100,
                  min.overlap = minsize,
                  extended.output = TRUE,
                  global.fdr = FALSE,
                  log.file = log)
