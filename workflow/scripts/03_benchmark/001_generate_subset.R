#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_files <- snakemake@input$scores
  output_file <- snakemake@output$output
}else{
  input_files <- list.files("results/03_benchmark/merged/01_input_bench", full.names = T)
}

## Libraries ---------------------------
library(tidyverse)

comb <- map_dfr(input_files, function(file){
  df <- read_csv(file, col_types = cols())
  df %>%
    pivot_longer(!experiment, values_to = "score", names_to = "kinase") %>%
    filter(!is.na(score)) %>%
    mutate(id = paste(experiment, kinase, sep = ":")) %>%
    dplyr::select(id)
})

shared <- comb %>%
  group_by(id) %>%
  summarise(n = n()) %>%
  filter(n == length(input_files)) %>%
  dplyr::select(id)

write_csv(shared, output_file)
