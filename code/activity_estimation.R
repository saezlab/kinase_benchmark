# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we estimate kinase
#' activities with all combinations of
#' KSN and method

library(tidyverse)
library(decoupleR)
library(tidyselect)
source("code/support_functions.R")
map(list.files("code", pattern = "run", full.names = T), source)

## Load data---------------------------
### phosphoproteomics
ccrcc <- read_csv(file = "data/CPTAC_estefani_clean_zip/ccrcc/ccrcc_phospho_tvalues.csv") %>%
  select(ID, NATvsTUM_t) %>%
  rename("ccrcc_tVal" = "NATvsTUM_t")
luad <- read_csv(file = "data/CPTAC_estefani_clean_zip/luad/luad_phospho_tvalues.csv") %>%
  select(ID, NATvsTUM_t) %>%
  rename("luad_tVal" = "NATvsTUM_t")
ucec <- read_csv(file = "data/CPTAC_estefani_clean_zip/ucec/ucec_phospho_tvalues.csv") %>%
  select(ID, NATvsTUM_t) %>%
  rename("ucec_tVal" = "NATvsTUM_t")

mat <- reduce(list(ccrcc, luad, ucec), full_join) %>%
  filter(rowSums(is.na(.)) == 0) %>%
  column_to_rownames("ID")


### networks
network_files <- list.files("data/networks", full.names = T)
networks <- map(network_files, read_csv)
names(networks) <-  str_remove(list.files("data/networks"), ".csv")

## ---------------------------
kin_scores <- map(networks, function(network){
  #prepare network
  network <- intersect_regulons(mat, network, .source = "source", .target = "target", minsize = 5)
  cor.source <- check_corr(network) %>% filter(correlation > 0.9) %>% pull(source.2)
  network <- network %>% filter(!source %in% cor.source)

  KARP <- run_KARP(mat, network)
  RoKAI_z <- run_zscore_RoKAI(mat, network)
  KSEA_z <- run_zscore_KSEA(mat, network)
  INKA <- run_INKA(mat, network)
  Rokai_lm <- run_lm_rokai(mat, network)
  decoupler <- decouple(mat = as.matrix(mat), network = network)

  decoupler <- decoupler %>%
    select(c(source, condition, score, statistic)) %>%
    rename("method" = "statistic")

  rbind(KARP, RoKAI_z, KSEA_z, INKA, Rokai_lm, decoupler)
})
