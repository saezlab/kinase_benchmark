if(exists("snakemake")){
  prior_file <- snakemake@input$prior
  method_file <- snakemake@input$method
  prior_out <- snakemake@output$prior_out
  method_out <- snakemake@output$method_out
}else{
  prior_file <- list.files("results/manuscript_figures/supp_files", full.names = T, pattern = "prior_comparison")
  method_file <- list.files("results/manuscript_figures/supp_files", full.names = T, pattern = "method_comparison")
  prior_out <- "results/manuscript_figures/supp_files/prior_comparison.csv"
  method_out <- "results/manuscript_figures/supp_files/method_comparison.csv"
}


## Libraries ---------------------------
library(tidyverse)

## Format data ---------------------------
prior_res <- map_dfr(prior_file, function(file_id){
    read_csv(file_id, col_types = cols()) %>% 
        add_column(benchmark = sub(".*_(.*?)\\.csv", "\\1", file_id), .before = 1)
    }) %>%
    mutate(benchmark = case_when(
        benchmark == "act" ~ "activating site-based",
        benchmark == "tumor" ~ "protein-based",
        benchmark == "perturbation" ~ "perturbation-based"
    )) %>%
    rename("evaluation approach" = "benchmark", "Library 1" = "comp1", "Library 2" = "comp2", "p-value" = "p.value", "adj. p-value" = "p.adj")


method_res <- map_dfr(method_file, function(file_id){
    read_csv(file_id, col_types = cols()) %>% 
        add_column(benchmark = sub(".*_(.*?)\\.csv", "\\1", file_id), .before = 1)
    }) %>%
    mutate(benchmark = case_when(
        benchmark == "act" ~ "activating site-based",
        benchmark == "tumor" ~ "protein-based",
        benchmark == "perturbation" ~ "perturbation-based"
    )) %>%
    rename("evaluation approach" = "benchmark", "Comp. method 1" = "comp1", "Comp. method 2" = "comp2", "p-value" = "p.value", "adj. p-value" = "p.adj")


write_csv(prior_res, prior_out)
write_csv(method_res, method_out)
