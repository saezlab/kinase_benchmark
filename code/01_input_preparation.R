# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we prepare the data for the different tools

library(tidyverse)
library(cmapR)
library(org.Hs.eg.db)

## Load data (containing 5 cell lines) ---------------------------
CPTAC_data <- list(ccrcc = read.csv("data/CPTAC_estefani_clean_zip/ccrcc/ccrcc_phospho_tvalues.csv"),
                   luad = read.csv("data/CPTAC_estefani_clean_zip/luad/luad_phospho_tvalues.csv"),
                   ucec = read.csv("data/CPTAC_estefani_clean_zip/ucec/ucec_phospho_tvalues.csv"))

CPTAC_data <- map(CPTAC_data, function(CPTAC_cohort){
  cbind(protein_symbol =  map_chr(str_split(CPTAC_cohort$ID, "_"), function(x) x[1]),
        pps =  map_chr(str_split(CPTAC_cohort$ID, "_"), function(x) x[2]),
        CPTAC_cohort)
})

# selection criteria for benchmark data?
benchmark_data <- readRDS("data/benchmark_data.rds")
benchmark_meta <- readRDS("data/benchmark_metaData.rds")

################# CPTAC ---------------------------
## Add uniprot ID, filter NAs in t-values, filter peptides with multipls pps
CPTAC_data_filtered <- map(CPTAC_data, function(CPTAC_cohort){
  mapped <- mapIds(org.Hs.eg.db,
                   keys = CPTAC_cohort$protein_symbol,
                   keytype = "SYMBOL",
                   column = c("UNIPROT"))

  CPTAC_cohort <- CPTAC_cohort %>%
    add_column(protein_uniprot = mapped, .after = "protein_symbol") %>%
    filter(rowSums(is.na(.)) == 0) %>%
    filter(str_count(pps, "[A-Z]") == 1)

  CPTAC_cohort %>% add_column(uniprot_ID = paste0(CPTAC_cohort$protein_uniprot, "_", CPTAC_cohort$pps), .after = "ID")

})

## RoKAI web application ---------------------------
map(names(CPTAC_data_filtered), function(cohort){
 cohort_data <- CPTAC_data_filtered[[cohort]]
  rokai_data <- data.frame(str_split_fixed(cohort_data$uniprot_ID, "_", 2)) %>%
    add_column(Quantification = cohort_data$NATvsTUM_logFC)

  colnames(rokai_data) <- c("Protein", "Position", "Quantification")

  write.csv(rokai_data, paste0("output/CPTAC/rokai/rokai_input_",cohort, ".csv"))
})


## KEA3 ---------------------------
# select proteins
map(names(CPTAC_data_filtered), function(cohort){
 cohort_data <- CPTAC_data_filtered[[cohort]]

  KEA3_data <-cohort_data %>%
    filter(NATvsTUM_adj.P.Val <= 0.05) %>%
    filter(abs(NATvsTUM_logFC) >= 2)

  KEA3_data <- unique(KEA3_data$protein_symbol)

  write.table(KEA3_data,  paste0("output/CPTAC/KEA3/KEA3_input_", cohort, ".txt"), row.names = F, col.names = F,  quote = FALSE)

})


## KSEA App ---------------------------
map(names(CPTAC_data_filtered), function(cohort){
 cohort_data <- CPTAC_data_filtered[[cohort]]
  KSEA_data <- data.frame(str_split_fixed(cohort_data$ID, "_", 2)) %>%
    add_column(p =cohort_data$NATvsTUM_adj.P.Val,
               FC = (2^cohort_data$NATvsTUM_logFC),
               Protein = str_split_fixed(cohort_data$uniprot_ID, "_", 2)[,1],
               Peptide = "NULL") %>%
    rename("Gene" = "X1", "Residue.Both" = "X2")

  KSEA_data <- KSEA_data[c(5,1,6,2,3,4)]

  write.csv(KSEA_data, paste0("output/CPTAC/KSEA/KSEA_input_",cohort, ".csv"), row.names = F)

})


## IKAP ---------------------------
# modified PSP dataset I downloaded for the KSEA App
map(names(CPTAC_data_filtered), function(cohort){
 cohort_data <- CPTAC_data_filtered[[cohort]]
  KSN_KSEA <- read.csv("data/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
  KSN_IKAP <- KSN_KSEA %>%
      filter(Source == "PhosphoSitePlus") %>%
      mutate(PSP = paste(SUB_GENE, SUB_MOD_RSD, sep = "_")) %>%
      select(c(GENE, PSP))


  IKAP_data <- data.frame(str_split_fixed(cohort_data$ID, "_", 2)[,1]) %>%
    add_column(peptide =cohort_data$ID,
      logFC =cohort_data$NATvsTUM_logFC)

  colnames(IKAP_data)[1] <- "protein"

  write.csv(KSN_IKAP, "code/IKAP/phosphositeplus.csv", row.names = F)
  write.csv(IKAP_data, paste0("output/CPTAC/IKAP/IKAP_input_CPTAC_",cohort, ".csv"), row.names = F)

})

## PTM-SEA ---------------------------
# create gct file
map(names(CPTAC_data_filtered), function(cohort){
 cohort_data <- CPTAC_data_filtered[[cohort]]
 cohort_data <- cohort_data %>%
    mutate(rid = paste0(str_replace(uniprot_ID, "_", ";"), "-p")) %>%
    select(c("rid", "NATvsTUM_t")) %>%
   rename("t" = "NATvsTUM_t") %>%
    column_to_rownames("rid")
 cohort_data_m <- as.matrix(cohort_data)

 cohort_data_gct <- new("GCT", mat=cohort_data_m)
 write_gct(cohort_data_gct, paste0("output/CPTAC/PTM-SEA/PTM-SEA_input_",cohort))
})
