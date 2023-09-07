#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_files <- snakemake@input$bench
  output_file <- snakemake@output$ov
  output_agg <- snakemake@output$ov_prior
}else{
  input_files <- list.files("results/hernandez/benchmark_files", pattern = "obs", full.names = T)
  output_file <- "results/hernandez/benchmark_res/overview/overview_bench.res"
  output_agg <- "results/hernandez/benchmark_res/overview/overview_bench_prior.res"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
overview_df <- map_dfr(input_files, function(file){
  df <- read_csv(file, col_types = cols())
  data.frame(n_experiments = nrow(df),
             n_kinases = str_split(df$perturb, ";") %>% unlist() %>% unique() %>% length(),
             method = str_remove(str_split(str_split(file, "/")[[1]][4], "-")[[1]][1], "obs_"),
             prior = str_remove(str_split(str_split(file, "/")[[1]][4], "-")[[1]][2], ".csv"))
}) %>%
  filter(!is.na(prior))

write_csv(overview_df, output_file)

write_csv(overview_df %>%
            group_by(prior) %>%
            summarise(kin = round(mean(n_kinases), digits = 1),
                      exp = round(mean(n_experiments), digits = 1)),
          output_agg)
