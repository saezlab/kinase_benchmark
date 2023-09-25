if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  GPS_file <- snakemake@input$GPS
  output_file_merged <- snakemake@output$merged
}else{
  ppsp_file <- "data/prior/phosphositeplus"
  GPS_file <- "data/prior/mmc4.xlsx"
  output_file_merged <- "results/GPSppsp.tsv"
}

## Libraries ---------------------------
library(tidyverse)
library(readxl)

## load table for kinase targets (GPS Gold Standard (GS) set) ---------------------------
GPS <- read_xlsx(GPS_file)

# filter to human kinases and extract 15mer centered on phosphosite; 266 kinases with at least 5 substrates
GPS <- GPS[GPS$Species=="Homo sapiens", ]

# Select 15mer
GPS$Site <- sapply(GPS$Sequence, function(x) paste0(unlist(strsplit(x, ""))[24:38], collapse = ""))
GPS$Site <- gsub("\\*","_",GPS$Site)

colnames(GPS)[1] <- "Kinase_UniProt"
colnames(GPS)[2] <- "Kinase"
colnames(GPS)[3] <- "Target_UniProt"

GPS <- GPS %>%
  select(Kinase, Kinase_UniProt, Target_UniProt, Position, Code, Site)

## Construct kinase-substrate interaction network ---------------------------
phosphositeplus <- read_tsv(ppsp_file, skip = 2)
phosphositeplus_human <- phosphositeplus %>%
  dplyr::filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>%
  mutate(Site = toupper(`SITE_+/-7_AA`))

phosphositeplus_human$Code <- stringr::str_extract(phosphositeplus_human$SUB_MOD_RSD, "^.{1}")
phosphositeplus_human$Position <- as.numeric(gsub("\\D", "", phosphositeplus_human$SUB_MOD_RSD))

phosphositeplus_human <- phosphositeplus_human %>%
  rename("Kinase" = GENE, "Kinase_UniProt" = KIN_ACC_ID, "Target_UniProt" = SUB_ACC_ID) %>%
  select(Kinase, Kinase_UniProt, Target_UniProt, Position, Code, Site) %>%
  mutate(Kinase_UniProt = map_chr(str_split(Kinase_UniProt, "-"), 1)) %>%
  mutate(Target_UniProt = map_chr(str_split(Target_UniProt, "-"), 1))


GPS_ppsp <- rbind(GPS, phosphositeplus_human) %>%
  distinct()

write_tsv(GPS_ppsp, output_file_merged)

