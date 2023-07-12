#' Merging decryptm metadata and matrix from different dataset

## Snakemake ---------------------------
if(exists("snakemake")){
  EC50_files <- snakemake@input$EC50_files
  meta_files <- snakemake@input$meta_files
  output_meta <- snakemake@output$output_meta
  output_EC50 <- snakemake@output$output_EC50
}else{
  EC50_files <- list.files("results/decryptm/phosphoproteome/", pattern = "^EC50", full.names = T)
  meta_files <- list.files("results/decryptm/phosphoproteome/", pattern = "^meta", full.names = T)
  output_meta <- "results/decryptm/processed_data/meta_data.csv"
  output_EC50 <- "results/decryptm/processed_data/EC50.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Merge metadata and EC50 matrices ---------------------------
meta <- map_dfr(meta_files, read_csv)
meta <- meta %>%
  mutate(sample = recode(sample,
                         "ddPTM_A549_Dasatinib_60min_R1_Rep_Dasatinib" = "ddPTM_A549_Dasatinib_60min_R1_Rep",
                         "ddPTM_A549_PD325901_60min_R1_Rep_PD325901" = "ddPTM_A549_PD325901_60min_R1_Rep",
                         "ddPTM_A549_Refametinib_60min_R1_Rep_Refametinib" = "ddPTM_A549_Refametinib_60min_R1_Rep",
                         "ddPTM_A549_Staursporin_60min_R1_Rep_Staursporin" = "ddPTM_A549_Staursporin_60min_R1_Rep",))

EC50_list <- map(EC50_files, read_csv)

EC50_df <- EC50_list %>% reduce(full_join, by = "pps_id")

colnames(EC50_df)[!colnames(EC50_df) %in% meta$sample]
meta$sample[!meta$sample %in% colnames(EC50_df)]

meta <- meta %>%
  filter(sample %in% colnames(EC50_df))

table(meta$drug, meta$cells)

## Save_data ---------------------------
write_csv(meta, output_meta)
write_csv(EC50_df, output_EC50)
