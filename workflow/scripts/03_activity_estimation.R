if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  dataset_name <- snakemake@wildcards$dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  output_file <- snakemake@output$rds
  scripts <- snakemake@input$scripts
  script_support <- snakemake@input$script_support
}else{
  dataset <- "data/CPTAC_phospho/luad_phospho_data_median_centered.tsv"
  dataset_name <- "luad"
  PKN <- "results/prior/iKiPdb.tsv"
  PKN_name <- "iKiPdb"
  output_file <- "results/activity_scores/luad_goldStandard.rds"
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
phospho <- read_table(dataset) %>%
    column_to_rownames("site")

## Prior knowledge Kinase-Substrate Networks
if (PKN_name %in% c("goldStandard", "combined", "NetworKIN")){
  net <- read.table(file = PKN, sep = "\t", header = T)
  colnames(net) <- c("source", "target")
  prior <- net %>%
    add_column(mor = 1)
} else {
  prior <- read.table(file = PKN, sep = "\t", header = T)
}

## Kinase activity estimation ---------------------------
results <- map_dfr(1:ncol(phospho), function(i){
  mat_i <- phospho[i] %>%
    drop_na()

  #prepare network
  prior_tmp <- intersect_regulons(mat_i, prior, .source = "source", .target = "target", minsize = 5)
  cor.source <- check_corr(prior_tmp) %>% filter(correlation > 0.9) %>% pull(source.2)

  if (dataset_name == "luad" & PKN_name == "iKiPdb"){
    prior_i <- prior %>% filter(!source == "EPHB3") %>% ungroup()
  } else {
    prior_i <- prior %>% filter(!source %in% cor.source) %>% ungroup()
  }


  # run activity estimation methods
  KARP <- run_KARP(mat_i, prior_i)
  RoKAI_z <- run_zscore_RoKAI(mat_i, prior_i)
  KSEA_z <- run_zscore_KSEA(mat_i, prior_i)
  INKA <- run_INKA(mat_i, prior_i)
  Rokai_lm <- run_lm_rokai(mat_i, prior_i)
  decoupler <- decouple(mat = as.matrix(mat_i), network = prior_i)
  fgsea <- run_fgsea(mat = as.matrix(mat_i), network = prior_i)
  #gsva <- run_gsva(mat = as.matrix(mat_i), network = prior_i)
  wmean <- run_wmean(mat = as.matrix(mat_i), network = prior_i)
  viper <- run_viper(mat = as.matrix(mat_i), network = prior_i)

  decoupler <- decoupler %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  fgsea <- fgsea %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  #gsva <- gsva %>%
  #  dplyr::select(c(source, condition, score, statistic)) %>%
  #  rename("method" = "statistic")
  wmean <- wmean %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")
  viper <- viper %>%
    dplyr::select(c(source, condition, score, statistic)) %>%
    dplyr::rename("method" = "statistic")

  results <- rbind(KARP, RoKAI_z, KSEA_z, INKA, Rokai_lm, decoupler, fgsea, #gsva,
                   wmean, viper)
})

results <- results %>%
  dplyr::filter(!method %in% c("wsum", "norm_wsum")) %>%
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
