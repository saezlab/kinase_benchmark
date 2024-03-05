#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input
  output_file <- snakemake@output
}else{
  input_file <- "results/hijazi/01_processed_data/pps_fifteenmer.csv"
  hernandez_file <- "results/hernandez/processed_data/pps_fifteenmer.csv"
  cptac_files <- list.files("data/CPTAC_phospho/final", full.names = T, pattern = "original")
  cptac_output <- "results/cptac/processed_data/pps_fifteenmer.csv"
  output_file <- "results/dataset/fifteenmer.csv"
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)
library(org.Hs.eg.db)

## Load phosphorylation sites ---------------------------
hijazi <- read_csv(input_file, col_types = cols())
hernandez <- read_csv(hernandez_file, col_types = cols())

## Prepare CPTAC ---------------------------
# Map targets to pps in data
pps <- map_dfr(cptac_files, function(file){
  df <- readRDS(file)
  data.frame(site = rownames(df))
})

pps <- pps %>%
  distinct(.keep_all = TRUE)

identifiers <- map_chr(str_split(pps$site,  "\\|"), 1)
identifiers <- map_chr(str_split(identifiers,  "\\."), 1)

pps_df <- data.frame(site = pps,
                     ensembl_gene_id = identifiers,
                     position = map_chr(str_split(pps$site,  "\\|"), 3),
                     surrounding = map_chr(str_split(pps$site,  "\\|"), 4))

## Convert Ensemble gene IDs to Gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name'),
             values = identifiers,
             mart = mart)

target_df <- full_join(pps_df, res, by = "ensembl_gene_id")

cptac <- data.frame(Protein = target_df$external_gene_name,
                    Position = gsub("[^0-9]", "", target_df$position),
                    Aminoacid = gsub("[0-9]", "", target_df$position),
                    Sequence = target_df$surrounding) %>%
  distinct()

write_csv(cptac, cptac_output)

## Merge and save dataset ---------------------------
full_df <- rbind(hijazi, hernandez, cptac) %>% distinct()
write_csv(full_df, output_file)
