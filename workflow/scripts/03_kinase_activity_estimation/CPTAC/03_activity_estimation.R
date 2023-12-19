if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  dataset_name <- snakemake@wildcards$dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  output_file <- snakemake@output$rds
  scripts <- snakemake@input$scripts
  script_support <- snakemake@input$script_support
}else{
  dataset <- "data/CPTAC_phospho/final/ccrcc_norm2prot_original_lm_log2_medCentRatio.rds"
  dataset_name <- "ccrcc"
  PKN <- "results/cptac/prior/ptmsigdb.tsv"
  PKN_name <- "ptmsigdb"
  output_file <- "results/cptac/activity_scores/ptmsigdb/original/original_brca-ptmsigdb.rds"
  scripts <- list.files("workflow/scripts/methods", pattern = "run", full.names = T)
  script_support <- "workflow/scripts/methods/support_functions.R"
}

## Libraries ---------------------------
library(decoupleR)
library(tidyselect)
library(tidyverse)
source(script_support)
map(scripts, source)

## Load data ---------------------------
### phosphoproteomics
phospho <- readRDS(dataset)

## Prior knowledge Kinase-Substrate Networks
prior <- read.table(file = PKN, sep = "\t", header = T)

## Kinase activity estimation ---------------------------
results <- map_dfr(1:ncol(phospho), function(i){
  mat_i <- phospho[,i] %>%
    as.data.frame()
  rownames(mat_i) <- rownames(phospho)
  mat_i <- mat_i %>%
    drop_na()
  colnames(mat_i) <- colnames(phospho)[i]

  #prepare network
  prior_tmp <- intersect_regulons(mat_i, prior, .source = "source", .target = "target", minsize = 5)
  cor.source <- check_corr(prior_tmp)
  filter_source <- cor.source %>% filter(correlation > 0.9) %>% pull(source.2) %>% unique()

  if (dataset_name == "luad" & PKN_name == "iKiPdb"){
    filter_source <- "EPHB3"
  }

  prior_i <- prior %>% filter(!source %in% filter_source) %>% ungroup()

  # run activity estimation methods
  KARP <- run_KARP(mat_i, prior)
  RoKAI_z <- run_zscore_RoKAI(mat_i, prior)
  KSEA_z <- run_zscore_KSEA(mat_i, prior)
  INKA <- run_INKA(mat_i, prior)
  Rokai_lm <- run_lm_rokai(mat_i, prior)
  mlm <- run_mlm(mat = as.matrix(mat_i), network = prior_i)
  ulm <- run_ulm(mat = as.matrix(mat_i), network = prior)
  wsum <- run_wsum(mat = as.matrix(mat_i), network = prior)
  fgsea <- run_fgsea(mat = as.matrix(mat_i), network = prior)
  wmean <- run_wmean(mat = as.matrix(mat_i), network = prior)
  viper <- run_viper(mat = as.matrix(mat_i), network = prior)

  # For mlm add correlated kinases again and assign the same score
  mlm <- mlm %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  recode_score <- map_dfr(filter_source, function(kin){
    cor.kinase <- cor.source %>%
      filter(source.2 == kin) %>%
      arrange(desc(correlation)) %>%
      pull(source)
    cor.kinase <- cor.kinase[1]

    score_kin <- mlm %>%
      filter(source == cor.kinase)
    score_kin$source <- kin
    score_kin
  })
  mlm <- rbind(mlm, recode_score) %>%
    distinct()

  # Rename decoupler output for rbind
  wsum <- wsum %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  ulm <- ulm %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  fgsea <- fgsea %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  wmean <- wmean %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  viper <- viper %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")

  results <- rbind(KARP, RoKAI_z, KSEA_z, INKA, Rokai_lm, ulm, mlm,
                   wsum, fgsea, wmean, viper) %>%
    filter(!is.infinite(score))
})

# Run methods implemented by Eric Kai from Zhang group
results_eric <- calculate_Kinase_Activity(mat = phospho, network = prior)
results_eric_long <- map_dfr(names(results_eric), function(method_i){
  tidyr::pivot_longer(results_eric[[method_i]] %>%
                        as.data.frame() %>%
                        rownames_to_column("source"),
                      !source,  names_to = "condition", values_to = "score") %>%
    add_column(method = method_i)
})

# Combine results
results <- rbind(results, results_eric_long)
results <- results %>%
  drop_na() %>%
  arrange(source) %>%
  distinct() %>%
  group_by(method)

results_list <- results %>%
  group_split()
names(results_list) <- group_keys(results)$method

res_wide <- map(results_list, function(result_method){
  result_method %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source")
})

saveRDS(res_wide, output_file)
