#'

## Snakemake ---------------------------
if(exists("snakemake")){
  phosphosite_file <- snakemake@input$phospho
  meta_file <- snakemake@input$meta
  bench_output <- snakemake@output$mat
  meta_output <- snakemake@output$meta_out
  exp_to_remove <- snakemake@params$msk_exp
}else{
  phosphosite_file <- "data/hernandez/annotations.xlsx"
  meta_file <- "data/hernandez/benchmark_data.xlsx"
  bench_output <- "results/hernandez/processed_data/benchmark_data.csv"
  meta_output <- "results/hernandez/processed_data/benchmark_metadata.csv"
  exp_to_remove <- c("705_225", "1298_272", "1288_272", "1291_272", "387_117", "1289_272", "1290_272", "1308_272", "699_225")
}

## Libraries ---------------------------
library(tidyverse)
library(readxl)

# Benchmark data
# load logFC and annotations
phosphosites <- read_excel(phosphosite_file, sheet = "Phosphosites")
logFC <- read_excel(phosphosite_file, sheet = "Foldchanges", na = "NA")

benchmark_data <- tibble(cbind(ID = paste0(phosphosites$gene_name, "|",
                                           phosphosites$residues, phosphosites$positions, "|",
                                           phosphosites$ensg, "|",
                                           phosphosites$ensp),
                               logFC))

# Filter out benchmamrk kinases not part of the KSN
metaData <- read_excel(meta_file, sheet = "KinaseConditionPairs")

# Meta data
metaData <- as_tibble(metaData) %>%
  dplyr::rename(id = Condition, target = Kinase, sign = Regulation) %>%
  filter(id %in% colnames(benchmark_data))
metaData <- metaData %>%
  mutate(sign = replace(sign, sign == "up", 1)) %>%
  mutate(sign = replace(sign, sign == "down", -1))
metaData$sign <- as.double(metaData$sign)

#manually change perturbation for 702_225 as this should be a knockin according to publication
metaData$sign[metaData$id == "702_225"] <- 1


## rename meta_data kinases
metaData <- metaData %>%
  mutate(target = recode(target,
                         "ABL" = "ABL1"))

metaData <- metaData %>%
  filter(!id %in% exp_to_remove)


# Save data
write_csv(benchmark_data, bench_output)
write_csv(metaData, meta_output)

