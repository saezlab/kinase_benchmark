# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will load the prior knowledge used in the benchmark


library(tidyverse)
library(OmnipathR)

## Omnipath---------------------------
### Construct kinase-substrate interaction network
omnipath_ptm <- get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation", "phosphorylation"), ]

# Filter out ProtMapper
omnipath_ptm_filtered <- omnipath_ptm %>%
  dplyr::filter(!(stringr::str_detect(omnipath_ptm$sources, "ProtMapper") & n_resources == 1))

# select target (substrate_genesymbol) and source (enzyme_genesymbol)
KSN <- omnipath_ptm_filtered[, c(4, 3, 7)]

# add phosphorylation site to target
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol, omnipath_ptm_filtered$residue_type, sep = "_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol, omnipath_ptm_filtered$residue_offset, sep = "")

KSN <- KSN %>% mutate(
  mor = case_when(
    modification == "phosphorylation" ~ 1,
    modification == "dephosphorylation" ~ -1
  )
)

# we remove duplicated edges
KSN$id <- paste(KSN$substrate_genesymbol, KSN$enzyme_genesymbol, sep = "")
KSN <- KSN[!(duplicated(KSN$id)), ]
KSN <- KSN[, -c(3,5)]

# rename KSN to fit decoupler format
names(KSN)[1:2] <- c("target", "source")

write_csv(KSN, "data/networks/omnipath.csv")


## PhosphoSitePlus---------------------------
### Construct kinase-substrate interaction network
omnipath_ptm <- get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation", "phosphorylation"), ]
omnipath_ptm_psp <- omnipath_ptm %>%
  dplyr::filter(stringr::str_detect(omnipath_ptm$sources, "PhosphoSite"))

# select target (substrate_genesymbol) and source (enzyme_genesymbol)
psp <- omnipath_ptm_psp[, c(4, 3, 7)]

# add phosphorylation site to target
psp$substrate_genesymbol <- paste(psp$substrate_genesymbol, omnipath_ptm_psp$residue_type, sep = "_")
psp$substrate_genesymbol <- paste(psp$substrate_genesymbol, omnipath_ptm_psp$residue_offset, sep = "")

psp <- psp %>% mutate(
  mor = case_when(
    modification == "phosphorylation" ~ 1,
    modification == "dephosphorylation" ~ -1
  )
)

# we remove duplicated edges
psp$id <- paste(psp$substrate_genesymbol, psp$enzyme_genesymbol, sep = "")
psp <- psp[!(duplicated(psp$id)), ]
psp <- psp[, -c(3,5)]

# rename KSN to fit decoupler format
names(psp)[1:2] <- c("target", "source")

write_csv(psp, "data/networks/phosphositeplus.csv")



## PTMsigDB ---------------------------
PTMsigDB <- gson::read.gmt("data/ptm.sig.db.all.uniprot.human.v1.9.0.gmt")
PTMsigDB.kinase <- PTMsigDB[str_detect(PTMsigDB$term, "KINASE"),]

df <- data.frame(source = map_chr(str_split(PTMsigDB.kinase$term, "_"), 2),
                 target = map_chr(str_split(PTMsigDB.kinase$gene, ";"), 1),
                 pps = str_remove(map_chr(str_split(PTMsigDB.kinase$gene, ";"), 2), "-p"))

df$target <-  map_chr(str_split(df$target, "-"), 1)

mapping <- read_tsv("data/uniprot_mapping.tsv") %>%
  rename("target" = From)

df <- df %>%
  left_join(mapping) %>%
  filter(!is.na(To))

df <- df %>%
  mutate(target = paste(df$To, df$pps, sep = "_")) %>%
  add_column(mor = 1) %>%
  select(source, target, mor)

write_csv(df, "data/networks/ptmsigdb.csv")

## NetworKIN ---------------------------
networkin <- read_csv("data/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv") %>%
  filter(Source == "NetworKIN")

df <- networkin %>%
  rename("source" = KINASE) %>%
  mutate(target = paste(networkin$SUB_GENE, networkin$SUB_MOD_RSD, sep = "_")) %>%
  add_column(mor = 1) %>%
  filter(networkin_score >= 5) %>%
  select(source, target, mor, networkin_score)

write_csv(df, "data/networks/networkin_5.csv")
