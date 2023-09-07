#' Option to scale the data between experiments

# Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  output_file <- snakemake@output$output
}else{
  input_file <- "results/hernandez/final_scores/ptmsigdb.rds"
  output_file <- "results/hernandez/final_scores/scaled/ptmsigdb.rds"
}

## Libraries ---------------------------
library(tidyverse)

method_tmp <- str_remove(str_split(output_file, "/")[[1]][4], ".csv")
input <- paste0("-",str_remove(str_split(input_file, "/")[[1]][4], ".rds"))
method <- gsub(input, "", method_tmp)

## Load scores and scale them ---------------------------
act_scores <- readRDS(input_file)

scaled_scores <- map(names(act_scores), function(mat_i){
  mat <- act_scores[[mat_i]]
  if (mat_i == "number_of_targets"){
    mat
  } else {
    scale(mat, center = FALSE, scale = TRUE)[,]
  }
})

names(scaled_scores) <- names(act_scores)

# Save scaled scores
saveRDS(scaled_scores, output_file)
