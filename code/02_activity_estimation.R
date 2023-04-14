# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we estimate kinase
#' activities with all combinations of
#' KSN and method

library(tidyverse)
library(decoupleR)
library(tidyselect)
map(list.files(file.path("code", "methods"), full.names = T), source)

## Load data ---------------------------
### phosphoproteomics
phospho.files <- list.files(path = "data/CPTAC_phospho", full.names = T)

phospho <- map(phospho.files, function(file){
  read_table(file) %>%
    column_to_rownames("site")})
names(phospho) <- str_remove(list.files(path = "data/CPTAC_phospho"), "_phospho_data_median_centered.tsv")

## Prior knowledge Kinase-Substrate Networks
network.files <- list.files(path = "data/prior", full.names = T)

networks <- map(network.files, function(file){
  net <- read.table(file = file, sep = "\t", header = T)
  colnames(net) <- c("source", "target")
  net %>%
    add_column(mor = 1)
})
names(networks) <- c("combined", "goldStandard", "NetworKin")


## Kinase activity estimation ---------------------------
scores <- map(names(networks), function(net_i){
  prior <- networks[[net_i]]
  kin_scores <- map(names(phospho), function(phospho_i){
    print(phospho_i)
    mat <- phospho[[phospho_i]]
    results <- map_dfr(1:ncol(mat), function(i){
      mat_i <- mat[i] %>%
        drop_na()

      #prepare network
      prior_i <- intersect_regulons(mat_i, prior, .source = "source", .target = "target", minsize = 5)
      cor.source <- check_corr(prior_i) %>% filter(correlation > 0.9) %>% pull(source.2)
      prior_i <- prior_i %>% filter(!source %in% cor.source) %>% ungroup()

      # run activity estimation methods
      KARP <- run_KARP(mat_i, prior_i)
      RoKAI_z <- run_zscore_RoKAI(mat_i, prior_i)
      KSEA_z <- run_zscore_KSEA(mat_i, prior_i)
      #INKA <- run_INKA(mat_i, prior_i)
      Rokai_lm <- run_lm_rokai(mat_i, prior_i)
      decoupler <- decouple(mat = as.matrix(mat_i), network = prior_i)
      fgsea <- run_fgsea(mat = as.matrix(mat_i), network = prior_i)
      #gsva <- run_gsva(mat = as.matrix(mat_i), network = prior_i)
      wmean <- run_wmean(mat = as.matrix(mat_i), network = prior_i)
      viper <- run_viper(mat = as.matrix(mat_i), network = prior_i)

      decoupler <- decoupler %>%
        select(c(source, condition, score, statistic)) %>%
        rename("method" = "statistic")
      fgsea <- fgsea %>%
        select(c(source, condition, score, statistic)) %>%
        rename("method" = "statistic")
      #gsva <- gsva %>%
      #  select(c(source, condition, score, statistic)) %>%
      #  rename("method" = "statistic")
      wmean <- wmean %>%
        select(c(source, condition, score, statistic)) %>%
        rename("method" = "statistic")
      viper <- viper %>%
        select(c(source, condition, score, statistic)) %>%
        rename("method" = "statistic")

      results <- rbind(KARP, RoKAI_z, KSEA_z, Rokai_lm, decoupler, fgsea, #gsva,
                       wmean, viper) %>%
        add_column(cancer = phospho_i)
    })

    results <- results %>%
      filter(!method %in% c("INKA_substrate_centric", "INKA_kinase_centric", "INKA", "wsum", "norm_wsum")) %>%
      group_by(method)

    results_list <- results %>%
      group_split()
    names(results_list) <- group_keys(results)$method

    res_wide <- map(results_list, function(result_method){
      result_method %>%
        select(source, condition, score) %>%
        pivot_wider(names_from = condition, values_from = score) %>%
        column_to_rownames("source")
    })
    saveRDS(res_wide, paste0("output/cptac/", phospho_i, "_", net_i,".rds"))
    res_wide
  })
  names(kin_scores) <- names(phospho)
  kin_scores
})

names(scores) <- names(networks)
