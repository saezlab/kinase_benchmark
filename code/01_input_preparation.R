# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we prepare the data for the different tools

library(tidyverse)
library(cmapR)

## Load data (containing 5 cell lines) ---------------------------
CRC_data <- readRDS("data/CRC_DE.rds")

## Add uniprot ID
met_crc_phosphoproteome <- as.data.frame(read_delim("../metformin_CRC/data/met_crc_phosphoproteome_raw.csv",
                                                    ",", escape_double = FALSE, trim_ws = TRUE))

IDs <- met_crc_phosphoproteome %>% select(ID, site.id.expanded.uniprot)

CRC_data <- map(CRC_data, function(cell_line_data){
  cell_line_data <- left_join(cell_line_data, IDs) %>% rename("Uniprot" = site.id.expanded.uniprot)
})

## RoKAI web application ---------------------------
map(names(CRC_data), function(cell_line){
  cell_line_data <- CRC_data[[cell_line]]
  rokai_data <- data.frame(str_split_fixed(cell_line_data$Uniprot, "_", 2)) %>%
    add_column(Quantification = cell_line_data$logFC)

  colnames(rokai_data) <- c("Protein", "Position", "Quantification")

  write.csv(rokai_data, paste0("output/rokai/rokai_input_", cell_line, ".csv"))
})


## KEA3 ---------------------------
# select proteins
map(names(CRC_data), function(cell_line){
  cell_line_data <- CRC_data[[cell_line]]

  KEA3_data <- cell_line_data %>%
    filter(adj.P.Val <= 0.05) %>%
    filter(abs(logFC) >= 2)

  KEA3_data <- data.frame(str_split_fixed(KEA3_data$ID, "_", 2))$X1

  write.table(KEA3_data,  paste0("output/KEA3/KEA3_input_", cell_line, ".txt"), row.names = F, col.names = F)

})


## KSEA App ---------------------------
map(names(CRC_data), function(cell_line){
  cell_line_data <- CRC_data[[cell_line]]
  KSEA_data <- data.frame(str_split_fixed(cell_line_data$ID, "_", 2)) %>%
    add_column(p = cell_line_data$adj.P.Val,
               FC = (2^cell_line_data$logFC),
               Protein = str_split_fixed(cell_line_data$Uniprot, "_", 2)[,1],
               Peptide = "NULL") %>%
    rename("Gene" = "X1", "Residue.Both" = "X2")

  KSEA_data <- KSEA_data[c(5,1,6,2,3,4)]

  write.csv(KSEA_data, paste0("output/KSEA/KSEA_input_", cell_line, ".csv"), row.names = F)

})


## IKAP ---------------------------
# modified PSP dataset I downloaded for the KSEA App
map(names(CRC_data), function(cell_line){
  cell_line_data <- CRC_data[[cell_line]]
  KSN_KSEA <- read.csv("data/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
  KSN_IKAP <- KSN_KSEA %>%
      filter(Source == "PhosphoSitePlus") %>%
      mutate(PSP = paste(SUB_GENE, SUB_MOD_RSD, sep = "_")) %>%
      select(c(GENE, PSP))


  IKAP_data <- data.frame(str_split_fixed(cell_line_data$ID, "_", 2)[,1]) %>%
    add_column(peptide = cell_line_data$ID,
      logFC = cell_line_data$logFC)

  colnames(IKAP_data)[1] <- "protein"

  write.csv(KSN_IKAP, "code/IKAP/phosphositeplus.csv", row.names = F)
  write.csv(IKAP_data, paste0("output/IKAP/IKAP_input_", cell_line, ".csv"), row.names = F)

})

## PTM-SEA ---------------------------
# create gct file
map(names(CRC_data), function(cell_line){
  cell_line_data <- CRC_data[[cell_line]]
  cell_line_data <- cell_line_data %>%
    mutate(rid = paste0(str_replace(Uniprot, "_", ";"), "-p")) %>%
    select(c("rid", "t")) %>%
    column_to_rownames("rid")
  cell_line_data_m <- as.matrix(cell_line_data)

  cell_line_data_gct <- new("GCT", mat=cell_line_data_m)
  write_gct(cell_line_data_gct, paste0("output/PTM-SEA/PTM-SEA_input_", cell_line))
})
