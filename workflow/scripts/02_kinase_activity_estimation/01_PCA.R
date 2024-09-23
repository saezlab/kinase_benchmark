if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  output_file <- snakemake@output$rds
  scripts_method <- snakemake@input$scripts
  remove_auto <- snakemake@params$rm_auto
  minsize <- snakemake@params$minsize
}else{
  dataset <- "results/01_processed_data/hernandez/data/benchmark_data.csv"
  PKN <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  PKN_name <- "phosphositeplus"
  output_file <- "results/02_activity_scores/hernandez/PCA/phosphositeplus.csv"
  remove_auto <- T
  scripts_method <- "workflow/scripts/methods/run_erics_methods.R"
  minsize <- 3
}

minsize <- as.double(minsize)

## Libraries ---------------------------
library(decoupleR)
library(tidyverse)
library(tidyselect)
source(scripts_method)

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

## Kinase activity estimation ---------------------------
results_wide <- calculate_Kinase_Activity(mat = phospho, network = prior, min_sites = minsize)
results_long <- map_dfr(names(results_wide), function(method_i){
  tidyr::pivot_longer(results_wide[[method_i]] %>%
                        as.data.frame() %>%
                        rownames_to_column("source"),
                      !source,  names_to = "condition", values_to = "score") %>%
    add_column(method = method_i)
}) %>%
  dplyr::filter(method == "PCA") %>%
  dplyr::filter(!is.na(score))

## Save results ---------------------------
write_csv(results_long, output_file)
