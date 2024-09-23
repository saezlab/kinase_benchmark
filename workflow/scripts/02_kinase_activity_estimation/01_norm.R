if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  output_file <- snakemake@output$rds
  remove_auto <- snakemake@params$rm_auto
  minsize <- snakemake@params$minsize
  cores <- snakemake@threads[[1]]
}else{
  dataset <- "results/01_processed_data/hernandez/data/benchmark_data.csv"
  PKN <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  PKN_name <- "phosphositeplus"
  output_file <- "results/02_activity_scores/hernandez/mean/ptmsigdb.csv"
  remove_auto <- T
  minsize <- 3
  cores <- 1
}

minsize <- as.double(minsize)

## Libraries ---------------------------
library(decoupleR)
library(tidyverse)
library(tidyselect)
library(furrr)

## Load data ---------------------------
### phosphoproteomics
phospho <- read_csv(dataset, col_types = cols()) %>%
  column_to_rownames("ID")

## Prior knowledge Kinase-Substrate Networks
prior <- read.table(file = PKN, sep = "\t", header = T) %>%
  dplyr::filter(str_detect(target, "\\|"))

remove_auto <- as.logical(remove_auto)
if (!remove_auto){
  prior <- prior %>%
    dplyr::mutate(target = str_remove(target, "\\|auto"))
}

#defining parallelisation
plan(multisession, workers = cores)

## Kinase activity estimation ---------------------------
results <- future_map_dfr(1:ncol(phospho), function(i){
  mat_i <- phospho[, i, drop = FALSE] %>%
    drop_na()

  # run activity estimation methods
  run_wmean(mat = as.matrix(mat_i), network = prior, minsize = minsize) %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
}) %>%
  dplyr::filter(method == "norm_wmean") %>%
  dplyr::mutate(method = recocde(method, norm_wmean = "norm_mean"))

## Save results ---------------------------
write_csv(results, output_file)
