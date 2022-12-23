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
