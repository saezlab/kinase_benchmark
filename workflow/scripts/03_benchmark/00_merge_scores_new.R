
## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  input_hernandez <- snakemake@input$hernandez
  input_tyrosine <- snakemake@input$tyr
  output_file <- snakemake@output$output
}else{
  input_file <- "results/02_activity_scores/hijazi/final_scores/GPS.rds"
  input_hernandez <- "results/02_activity_scores/hernandez/final_scores/GPS.rds"
  input_tyrosine <- "results/02_activity_scores/tyrosine/final_scores/GPS.rds"
  output_file <- "results/02_activity_scores/hernandez/merged/GPS.rds"
}

## Libraries ---------------------------
library(tidyverse)

## Load scores and meta ---------------------------
act_scores_hijazi <- readRDS(input_file)
act_scores_hernandez <- readRDS(input_hernandez)
act_scores_tyrosine <- readRDS(input_tyrosine)

## merge datasets
act_scores <- map(names(act_scores_hijazi), function(method_idx){
  hijazi <- act_scores_hijazi[[method_idx]] %>%
    as.data.frame() %>%
    rownames_to_column("kinase")
  hernandez <- act_scores_hernandez[[method_idx]] %>%
    as.data.frame() %>%
    rownames_to_column("kinase")
  tyrosine <- act_scores_tyrosine[[method_idx]] %>%
    as.data.frame() %>%
    rownames_to_column("kinase")

  full_join(hijazi, hernandez, by = "kinase") %>%
    full_join(tyrosine, by = "kinase") %>%
    column_to_rownames("kinase")
})
names(act_scores) <- names(act_scores_hijazi)

## Save merged data ---------------------------
saveRDS(act_scores, output_file)
