if(exists("snakemake")){
  GPS_file <- snakemake@input$GPS_file
  file_datasets <- snakemake@input$file_dataset
  decryptm_dataset <- snakemake@input$decryptm
  output_file <- snakemake@output$tsv
  output_file_decryptm <- snakemake@output$out_decryptm
}else{
  GPS_file <- "data/prior/mmc4.xlsx"
  output_file <- "results/prior/GPS.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
  decryptm_dataset <- "results/decryptm/processed_data/R2_pEC50.csv"
  output_file_decryptm <- "results/decryptm/prior/GPS.tsv"
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

#site info to use for mapping to rowname:
GPS$Site1 <- paste0(GPS$Site, "|", GPS$`UniProt ID...3`)

## Mapping prior to rownames in data ---------------------------
#map this id to Uniprot/Swissprot id used for GPS GS set; first establish ensembl > Swissprot key
mappings <- read.table("data/knowledge_based_isoform_selection_v1.4.txt", sep = "\t", stringsAsFactors = F, header = T)
#mappings <- mappings[!duplicated(paste0(mappings$gene, mappings$protein, mappings$SwissProt)), ]
mappings_sp <- mappings[!is.na(mappings$SwissProt), ]
mappings_sp <- mappings_sp[!duplicated(paste0(mappings_sp$protein, mappings_sp$SwissProt)), ]
mappings_sp$protein <- sub("\\..*", "", mappings_sp$protein)
rownames(mappings_sp) <- mappings_sp$protein

## Prepare CPTAC ---------------------------
# Prepare phosphorylation sites in data sets
pps <- map_dfr(file_datasets, function(file){
  df <- read_tsv(file, col_types = cols())
  data.frame(site = df$site)
})

pps <- pps %>%
  distinct(.keep_all = TRUE)

pps$fifteenmer <- sub("^[^|]*\\|[^|]*\\|[^|]*\\|", "", pps$site)
pps$fifteenmer <- sub("\\|.*", "", pps$fifteenmer)
pps$prot <- sub("^[^|]*\\|", "", pps$site)
pps$prot <- gsub("\\|.*", "", pps$prot)
pps$prot <- gsub("\\..*", "", pps$prot)
sites_sp <- pps$prot[pps$prot %in% mappings_sp$protein]
pps$sp <- NA
pps[pps$prot %in% sites_sp, "sp"] <- mappings_sp[pps[pps$prot %in% sites_sp, "prot"], "SwissProt"]
pps$Site1 <- paste0(pps$fifteenmer, "|", pps$sp)

## Prepare decryptm ---------------------------
decryptm_df <- read_csv(decryptm_dataset)
decryptm_identifiers <- decryptm_df %>%
  dplyr::select(pps_id) %>%
  mutate(external_gene_name = map_chr(str_split(pps_id, "\\|"), 2))  %>%
  mutate(position = map_chr(str_split(pps_id, "\\|"), 4)) %>%
  dplyr::rename("site" = pps_id)

## Convert Ensemble gene IDs to Gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('external_gene_name',
                            'uniprot_gn_id'),
             values = decryptm_identifiers$external_gene_name,
             mart = mart)

decryptm_ids <- full_join(decryptm_identifiers, res, by = "external_gene_name", relationship = "many-to-many") %>%
  mutate(Site1 = paste0(position, "|", uniprot_gn_id))

## merge with network ---------------------------
## CPTAC ---------------------------
GPS_mapped <- left_join(GPS, pps, by = "Site1", relationship = "many-to-many") %>%
  mutate(target = case_when(
    is.na(site) ~ Site1,
    !is.na(site) ~ site
  )) %>%
  mutate(target = case_when(
    `UniProt ID...1` == `UniProt ID...3` ~ paste0(target, "|auto"), #mark autophosphorylation
    `UniProt ID...1` != `UniProt ID...3` ~ target
  )) %>%
  dplyr::rename("source" = Kinase) %>%
  dplyr::select(source, target) %>%
  add_column(mor = 1) %>%
  distinct()

## Save prior
write_tsv(GPS_mapped, output_file)

## decryptm ---------------------------
GPS_mapped <- left_join(GPS, decryptm_ids, by = "Site1", relationship = "many-to-many") %>%
  mutate(target = case_when(
    is.na(site) ~ Site1,
    !is.na(site) ~ site
  )) %>%
  mutate(target = case_when(
    `UniProt ID...1` == `UniProt ID...3` ~ paste0(target, "|auto"), #mark autophosphorylation
    `UniProt ID...1` != `UniProt ID...3` ~ target
  )) %>%
  dplyr::rename("source" = Kinase) %>%
  dplyr::select(source, target) %>%
  add_column(mor = 1) %>%
  distinct()

## Save prior
write_tsv(GPS_mapped, output_file_decryptm)
