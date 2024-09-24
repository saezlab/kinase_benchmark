if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  output_file <- snakemake@output$rds
  scripts_method <- snakemake@input$scripts
  script_support <- snakemake@input$script_support
  remove_auto <- snakemake@params$rm_auto
  minsize <- snakemake@params$minsize
  background <- snakemake@params$background
  cores <- snakemake@threads[[1]]
}else{
  dataset <- "results/01_processed_data/hernandez/data/benchmark_data.csv"
  PKN <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  PKN_name <- "phosphositeplus"
  output_file <- "results/02_activity_scores/hernandez/chisq/ptmsigdb.csv"
  scripts_method <- "workflow/scripts/methods/run_chisq.R"
  script_support <- "workflow/scripts/methods/support_functions.R"
  remove_auto <- T
  minsize <- 3
  cores <- 1
  background <- 20000
}

minsize <- as.double(minsize)

## Libraries ---------------------------
library(decoupleR)
library(tidyverse)
library(tidyselect)
library(furrr)
source(scripts_method)
source(script_support)

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

  n_up <- sum(mat_i[ , 1] >= 1)
  n_bottom <- sum(mat_i[ , 1] <= -1)

  # run activity estimation methods
  up <- run_chisq(mat = as.matrix(mat_i), network = prior, minsize = minsize, n_background = background, n_up = n_up, n_bottom = 0) %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")

  down <- run_chisq(mat = as.matrix(mat_i), network = prior, minsize = minsize, n_background = background, n_bottom = n_bottom, n_up = 0) %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic") %>%
    dplyr::mutate(score = -score)

  df <- rbind(up, down) %>%
    group_by(source, condition) %>%
    dplyr::filter(abs(score) == max(abs(score))) %>%# choose activity score with higher value
    ungroup() %>%
    distinct()
  
  df[!duplicated(df[, c("source", "condition")]),]

})

## Save results ---------------------------
write_csv(results, output_file)
