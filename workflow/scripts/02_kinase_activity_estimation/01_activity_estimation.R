if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  output_file <- snakemake@output$rds
  scripts <- snakemake@input$scripts
  script_support <- snakemake@input$script_support
  remove_auto <- snakemake@params$rm_auto
  minsize <- snakemake@params$minsize
  cores <- snakemake@threads[[1]]
}else{
  dataset <- "results/01_processed_data/hernandez/data/benchmark_data.csv"
  PKN <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  PKN_name <- "phosphositeplus"
  output_file <- "results/hijazi/03_activity_scores/ptmsigdb.rds"
  remove_auto <- T
  scripts <- list.files("workflow/scripts/methods", pattern = "run", full.names = T)
  script_support <- "workflow/scripts/methods/support_functions.R"
  minsize <- 3
  cores <- 1
}

minsize <- as.double(minsize)

## Libraries ---------------------------
library(decoupleR)
library(tidyselect)
library(tidyverse)
library(furrr)
source(script_support)
map(scripts, source)

## Load data ---------------------------
### phosphoproteomics
phospho <- read_csv(dataset, col_types = cols()) %>%
  column_to_rownames("ID")

## Prior knowledge Kinase-Substrate Networks
prior <- read.table(file = PKN, sep = "\t", header = T)

remove_auto <- as.logical(remove_auto)
if (!remove_auto){
  prior <- prior %>%
    dplyr::mutate(target = str_remove(target, "\\|auto"))
}

#defining parallelisation
plan(multisession, workers = cores)

## Kinase activity estimation ---------------------------
results <- future_map_dfr(1:ncol(phospho), function(i){
  print(paste0(i, "/", ncol(phospho)))
  if(str_detect(dataset, "original")){
    mat_i <- phospho[i] %>%
      drop_na()
  } else {
    mat_i <- phospho[,i] %>%
      data.frame() %>%
      drop_na()

    colnames(mat_i) <- colnames(phospho)[i]
  }

  mat_i <- phospho[i] %>%
    drop_na()

  # run activity estimation methods
  KARP <- run_KARP(mat_i, prior, minsize = minsize)
  RoKAI_z <- run_zscore_RoKAI(mat_i, prior, minsize = minsize)
  KSEA_z <- run_zscore_KSEA(mat_i, prior, minsize = minsize)
  INKA <- run_INKA(mat_i, prior, kinase_mapping = F, minsize = minsize)
  Rokai_lm <- run_lm_rokai(mat_i, prior, minsize = minsize)
  ulm <- run_ulm(mat = as.matrix(mat_i), network = prior, minsize = minsize)
  wsum <- run_wsum(mat = as.matrix(mat_i), network = prior, minsize = minsize)
  fgsea <- run_fgsea(mat = as.matrix(mat_i), network = prior, minsize = minsize)
  wmean <- run_wmean(mat = as.matrix(mat_i), network = prior, minsize = minsize)
  viper <- run_viper(mat = as.matrix(mat_i), network = prior, minsize = minsize)

  run.mlm.cor <- function(mat, network) {
    tmp <- FALSE
    #prepare network
    prior_tmp <- intersect_regulons(mat, network, .source = "source", .target = "target", minsize = minsize)
    cor.source <- check_corr(prior_tmp)
    i <- 0

    while (any(class(tmp) == "error") | any(class(tmp) == "logical")) {
      if (i == 0){
        remove <- 0.9
      } else {
        remove <- unique(cor.source$correlation)[i]
      }

      filter_source <- cor.source %>% filter(correlation >= remove) %>% pull(source.2) %>% unique()

      prior_i <- network %>% filter(!source %in% filter_source) %>% ungroup()

      tmp <- tryCatch(run_mlm(mat = as.matrix(mat), network = prior_i, minsize = minsize),
                      error=function(e) e, warning=function(w) w)

      i <- i + 1
    }

    mlm <- run_mlm(mat = as.matrix(mat_i), network = prior_i, minsize = minsize)
    return(list(score = mlm, filter_kin = filter_source, cor = cor.source))
  } # test correlation sequentially
  mlm_list <- run.mlm.cor(mat = mat_i, network = prior)
  mlm <- mlm_list$score
  filter_source <- mlm_list$filter_kin
  cor.source <- mlm_list$cor
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

 rbind(KARP, RoKAI_z, KSEA_z, INKA, Rokai_lm, ulm, mlm,
                   wsum, fgsea, wmean, viper)
})

# Run methods implemented by Eric Kai from Zhang group
results_eric <- calculate_Kinase_Activity(mat = phospho, network = prior, min_sites = minsize)
results_eric_long <- future_map_dfr(names(results_eric), function(method_i){
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
