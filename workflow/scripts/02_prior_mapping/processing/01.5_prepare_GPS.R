if(exists("snakemake")){
  GPS_file <- snakemake@input$GPS_file
  output_file <- snakemake@output$tsv
}else{
  GPS_file <- "data/prior/mmc4.xlsx"
  output_file <- "results/prior/raw/GPS.tsv"
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
GPS_mapped <- GPS %>%
  dplyr::select(Kinase, target_site, Site) %>%
  dplyr::rename("source" = Kinase, "target" = target_site, "sequence" = Site) %>%
  mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
  mutate(position = map_chr(str_split(target, "_"), 2)) %>%
  mutate(mor = 1) %>%
  mutate(source = recode(source,
                         "ADRBK1" = "GRK2",
                         "GSG2" = "HASPIN",
                         "ICK" = "CILK1",
                         "MLTK" = "MAP3K20",
                         "MST4" = "STK26",
                         "NIM1" = "NIM1K",
                         "PAK7" = "PAK5",
                         "MPS1" = "TTK"
  )) %>%
  dplyr::select(source, target, target_protein, position, mor, sequence) %>%
  distinct()

## Save prior
write_tsv(GPS_mapped, output_file)
