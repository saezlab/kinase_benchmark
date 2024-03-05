print("firsttest")
if(exists("snakemake")){
  scores_ptmsea_file <- snakemake@input$file_ptmsea
  scores_file <- snakemake@input$file_scores
  dataset_name <- snakemake@wildcards$dataset
  PKN_name <- snakemake@wildcards$PKN
  methods_to_remove <- snakemake@params$rm_methods
  output <- snakemake@output$rds
}else{
  scores_ptmsea_file <- "results/activity_scores_ptmsea/hnscc_GPS-scores.gct"
  scores_file <- "results/activity_scores/GPS/hnscc_GPS.rds"
  dataset_name <- "hnscc"
  PKN_name <- "GPS"
  output <- "results/final_scores/GPS/hnscc_GPS.rds"
  methods_to_remove <- c("corr_wmean", "corr_wsum", "norm_wsum", "wmean")
}


## Libraries ---------------------------
library(tidyverse)

## load scores ---------------------------
ptmsea <- read.delim(file=scores_ptmsea_file, skip=2)
ptmsea <- ptmsea %>%
  dplyr::select(colnames(ptmsea)[!str_detect(colnames(ptmsea), "Signature")]) %>%
  dplyr::select(-No.columns.scored) %>%
  column_to_rownames("id")


## combine scores ---------------------------
scores <- readRDS(scores_file)
scores$ptmsea <- ptmsea

scores <- scores[!names(scores) %in% methods_to_remove]

## save final scores ---------------------------
saveRDS(scores, file = output)

