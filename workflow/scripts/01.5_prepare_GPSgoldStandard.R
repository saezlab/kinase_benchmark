if(exists("snakemake")){
  GPS_file <- snakemake@input$GPS_file
  output_file <- snakemake@output$tsv
  file_datasets <- snakemake@input$file_dataset
}else{
  GPS_file <- "data/prior/mmc4.xlsx"
  output_file <- "results/prior/GPS_goldStandard.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
}

## Libraries ---------------------------
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


# Mapping
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
  add_column(mor = 1)

## Save prior ---------------------------
write_tsv(GPS_mapped, output_file)
