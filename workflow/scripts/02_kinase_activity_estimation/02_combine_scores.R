if(exists("snakemake")){
  score_files <- snakemake@input$file_scores
  output <- snakemake@output$rds
}else{
  score_files <- list.files("results/02_activity_scores/hernandez", recursive = T, full.names = T, pattern = "phosphositeplus.csv")
  output <- "results/hernandez/scores/GPS.rds"
}

## Libraries ---------------------------
library(tidyverse)

## load scores ---------------------------
results <- map_dfr(score_files, ~ read_csv(.x, col_types = cols())) %>%
  dplyr::select(source, condition, score, method)

## combine scores ---------------------------
results <- results %>%
  drop_na() %>%
  arrange(source) %>%
  distinct() %>%
  group_by(method)

results_list <- results %>%
  group_split()
names(results_list) <- group_keys(results)$method

res_wide <- map(results_list, function(result_method){
  result_method %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source")
})

## save final scores ---------------------------
saveRDS(res_wide, output)

