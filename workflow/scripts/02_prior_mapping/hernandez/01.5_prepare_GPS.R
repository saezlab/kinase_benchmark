if(exists("snakemake")){
  GPS_file <- snakemake@input$GPS_file
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  GPS_file <- "data/prior/mmc4.xlsx"
  output_file <- "results/hernandez/prior/GPS.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
}

## Libraries ---------------------------
library(biomaRt)
library(readxl)
library(tidyverse)

## load table for kinase targets (GPS Gold Standard (GS) set) ---------------------------
GPS <- read_xlsx(GPS_file)

# filter to human kinases and extract 15mer centered on phosphosite; 266 kinases with at least 5 substrates
GPS <- GPS[GPS$Species=="Homo sapiens", ]

# Select 15mer
GPS$Site <- sapply(GPS$Sequence, function(x) paste0(unlist(strsplit(x, ""))[24:38], collapse = ""))
GPS$Site <- gsub("\\*","_",GPS$Site)

colnames(GPS)[2] <- "Kinase"
colnames(GPS)[3] <- "Target"

GPS <- GPS %>%
  mutate(Target = map_chr(str_split(Target, "-"), 1))

# Transfer target
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res_kin <- getBM(attributes = c('uniprot_gn_id',
                                'external_gene_name'),
                 values = GPS$Target,
                 mart = mart)  %>%
  dplyr::rename("Target" = uniprot_gn_id) %>%
  dplyr::rename("target_gene_name" = external_gene_name)

GPS <- left_join(GPS, res_kin, by = "Target", relationship = "many-to-many") %>%
  filter(!is.na(target_gene_name))
GPS <- GPS %>%
  mutate(target_site = paste0(target_gene_name, "_", Code, Position))

## Prepare data ---------------------------
hernandez_df <- read_csv(dataset)
identifiers <- hernandez_df %>%
  dplyr::select(ID) %>%
  dplyr::mutate(gene_name = map_chr(str_split(ID, "\\|"), 1))  %>%
  dplyr::mutate(position = map_chr(str_split(ID, "\\|"), 2)) %>%
  dplyr::mutate(ensg = map_chr(str_split(ID, "\\|"), 3)) %>%
  dplyr::mutate(ensp = map_chr(str_split(ID, "\\|"), 4)) %>%
  dplyr::mutate(target_site = paste0(gene_name, "_", position)) %>%
  dplyr::rename("site" = ID)

## merge with network ---------------------------
GPS_mapped <- left_join(GPS, identifiers, by = "target_site", relationship = "many-to-many") %>%
  mutate(target = case_when(
    is.na(site) ~ target_site,
    !is.na(site) ~ site
  )) %>%
  mutate(target = case_when(
    str_detect(pattern = Kinase, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = Kinase, string = target_site) ~ target
  )) %>%
  dplyr::rename("source" = Kinase) %>%
  dplyr::select(source, target) %>%
  add_column(mor = 1) %>%
  distinct()

## Save prior
write_tsv(GPS_mapped, output_file)
