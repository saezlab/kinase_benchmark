if(exists("snakemake")){
  networkin_file <- snakemake@input$networkin_file
  output_file <- snakemake@output$tsv
  networkin_score <- snakemake@params$score
}else{
  networkin_file <- "data/prior/networkin_human_predictions_3.1.tsv"
  output_file <- "results/prior/raw/networkin.tsv"
  networkin_score <- 5
}
networkin_score <- as.numeric(networkin_score)

## Libraries ---------------------------
library(tidyverse)

## load NetworKIN ---------------------------
# load NetworKin prediction scores; convert kinase names to gene symbol and filter to kinases in KS table
nwkin <- read.table(networkin_file, sep = "\t",  stringsAsFactors = F, header= T, comment.char = "")

# rename kinases
nwkin$id <- toupper(nwkin$id)
nwkin <- nwkin %>%
  mutate(id = recode(id,
                     "ABL" = "ABL1",
                     "AURORAA" = "AURKA",
                     "AURORAA" = "AURKA",
                     "AURORAB" = "AURKB",
                     "AURORAC" = "AURKC",
                     "BRK" = "PTK6",
                     "CAMKIALPHA" = "CAMK1",
                     "CAMKIIALPHA" = "CAMK2A",
                     "CAMKIIBETA" = "CAMK2B",
                     "CAMKIIDELTA" = "CAMK2D",
                     "CAMKIIGAMMA" = "CAMK2G",
                     "CAMKIV" = "CAMK4",
                     "CAMKIDELTA" = "CAMK1D",
                     "MRCKB" = "CDC42BPB",
                     "DMPK2" = "CDC42BPG",
                     "ICK" = "CILK1",
                     "CK1ALPHA" = "CSNK1A1",
                     "CK1DELTA" = "CSNK1D",
                     "CK1EPSILON" = "CSNK1E",
                     "CK1GAMMA1" = "CSNK1G1",
                     "CK1GAMMA2" = "CSNK1G2",
                     "CK1GAMMA3" = "CSNK1G3",
                     "CK2ALPHA" = "CSNK2A1",
                     "CK2A2" = "CSNK2A2",
                     "DNAPK" = "PRKDC",
                     "GSK3ALPHA" = "GSK3A",
                     "GSK3BETA" = "GSK3B",
                     "IKKALPHA" = "CHUK",
                     "IKKBETA" = "IKBKB",
                     "IRR" = "INSRR",
                     "LKB1" = "STK11",
                     "MAP4K6" = "MINK1",
                     "MRCKA" = "CDC42BPA",
                     "RON" = "MST1R",
                     "MST2" = "STK3",
                     "MRCKA" = "CDC42BPA",
                     "TRKC" = "NTRK3",
                     "PDGFRBETA" = "PDGFRB",
                     "PDGFRALPHA" = "PDGFRA",
                     "PDHK1" = "PDK1",
                     "PDHK2" = "PDK2",
                     "PDHK3" = "PDK3",
                     "PDHK4" = "PDK4",
                     "PKBALPHA" = "AKT1",
                     "PKBBETA" = "AKT2",
                     "PKBGAMMA" = "AKT3",
                     "PKCALPHA" = "PRKCA",
                     "PKCEPSILON" = "PRKCE",
                     "PKCETA" = "PRKCH",
                     "PKCGAMMA" = "PRKCG",
                     "PKCIOTA" = "PRKCI",
                     "PKCTHETA" = "PRKCQ",
                     "PKG1CGKI" = "PRKG1",
                     "AMPKA1" = "PRKAA1",
                     "AMPKA2" = "PRKAA2",
                     "RSK3" = "RPS6KA2",
                     "P70S6K" = "RPS6KB2",
                     "LOK" = "STK10",
                     "MST3" = "STK24",
                     "YSK1" = "STK25",
                     "MST4" = "STK26",
                     "TRKA" = "NTRK1",
                     "TRKB" = "NTRK2",
                     "CHAK2" = "TRPM6",
                     "YES" = "YES1",
                     "PKAALPHA" = "PRKACA",
                     "PKABETA" = "PRKACB",
                     "PKAGAMMA" = "PRKACG"))

# Select targets above certain NetworKIN score
nwkin_filtered <- nwkin[nwkin$networkin_score >= networkin_score, ]
nwkin_filtered$AA <- unlist(str_extract_all(nwkin_filtered$sequence, '[[:lower:]]+'))

nwkin_filtered <- nwkin_filtered %>%
  mutate(target = map_chr(str_split(X.substrate, " "), 1)) %>%
  mutate(target_site = paste0(target, "_", toupper(AA), position))


## Merge with network ---------------------------
nwkin_df <- nwkin_filtered %>%
  dplyr::select(id, sequence, target_site) %>%
  dplyr::rename("source" = id, "target" = target_site) %>%
  mutate(sequence = toupper(sequence)) %>%
  mutate(sequence = str_replace_all(sequence, "-", "_")) %>%
  mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
  mutate(position = map_chr(str_split(target, "_"), 2)) %>%
  mutate(mor = 1) %>%
  dplyr::select(source, target, target_protein, position, mor, sequence) %>%
  distinct()

## Save prior
write_tsv(nwkin_df, output_file)


