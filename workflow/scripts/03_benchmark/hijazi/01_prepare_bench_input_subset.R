#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  input_hernandez <- snakemake@input$hernandez
  meta_file <- snakemake@input$meta
  meta_hernandez <- snakemake@input$metah
  output_file <- snakemake@output$output
  meta_out <- snakemake@output$meta_out
  rm_experiments <- snakemake@params$rm_exp
}else{
  input_file <- list.files("results/hijazi/04_final_scores/scaled", full.names = T, pattern = "phosphositeplus")
  input_hernandez <- list.files("results/hernandez/final_scores/scaled", full.names = T, pattern = "phosphositeplus")
  meta_file <- "results/hijazi/01_processed_data/benchmark_metadataPrior.csv"
  meta_hernandez <- "results/hernandez/processed_data/benchmark_metadata.csv"
  methods <- c("KARP", "KS", "KSEA_z", "PC1", "RoKAI_z", "UQ", "Wilcox", "fgsea", "mean", "median", "mlm", "norm_fgsea","norm_wmean","number_of_targets","rokai_lm","ulm","viper","wsum", "ptmsea")
  rm_experiments <- F
}

## Libraries ---------------------------
library(tidyverse)

## Load data ---------------------------
priors <- map_chr(str_split(input_file, "/"),5) %>% str_remove(".rds")
priors <- priors[map_lgl(priors, function(x) (any(str_detect(input_file, x))) & any(str_detect(input_hernandez, x)))]

act_scores_list <- map(priors, function(prior_idx){
  act_scores_hijazi <- input_file[str_detect(input_file, paste0("/", prior_idx, ".rds"))] %>% readRDS()
  act_scores_hernandez <- input_hernandez[str_detect(input_hernandez, paste0("/", prior_idx, ".rds"))] %>% readRDS()

  ## merge datasets
  act_scores <- map(names(act_scores_hijazi), function(method_idx){
    hijazi <- act_scores_hijazi[[method_idx]] %>%
      as.data.frame() %>%
      rownames_to_column("kinase")
    hernandez <- act_scores_hernandez[[method_idx]] %>%
      as.data.frame() %>%
      rownames_to_column("kinase")

    full_join(hijazi, hernandez, by = "kinase") %>%
      column_to_rownames("kinase")
  })
  names(act_scores) <- names(act_scores_hijazi)
  act_scores
})

names(act_scores_list) <- priors

## Prepare meta data
obs_hijazi <- read_csv(meta_file, col_types = cols()) %>%
  dplyr::select(id, sign, target)
obs_hernandez <- read_csv(meta_hernandez, col_types = cols()) %>%
  dplyr::select(id, sign, target)
obs <- rbind(obs_hijazi, obs_hernandez)

# Set perturb
obs_targets <- obs %>%
  dplyr::rename("perturb" = target)

# filter out experiments with unknown target (e.g. several members)
target_df <- obs_targets %>%
  group_by(id) %>%
  summarise(targets = paste(perturb, collapse = ";"), sign = unique(sign)) %>%
  dplyr::rename("perturb" = targets) %>%
  dplyr::filter(!perturb == "" | !is.na(perturb)) %>%
  dplyr::rename("sample" = id)

## Get subset ---------------------------
scores_long <- map_dfr(names(act_scores_list), function(score_i){
  map_dfr(methods, function(method_i){
    act_scores_list[[score_i]][[method_i]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      add_column(prior = score_i) %>%
      add_column(method = method_i) %>%
      filter(!is.na(score))
  })
})

th_all <- length(unique(scores_long$method)) * length(unique(scores_long$prior))

msk <- scores_long %>%
  group_by(kinase, sample) %>%
  summarise(n = n()) %>%
  filter(n >= th_all) %>%
  mutate(tmp = paste(kinase, sample, sep = ":")) %>%
  pull(tmp)

scores_filtered <- scores_long %>%
  mutate(tmp = paste(kinase, sample, sep = ":")) %>%
  filter(tmp %in% msk) %>%
  dplyr::select(-tmp)

## Prepare input ---------------------------
map(priors, function(prior_i){
  act_scores <- map(methods, function(method_idx){
    scores_filtered %>%
      dplyr::filter(prior == prior_i) %>%
      filter(method == method_idx) %>%
      dplyr::select(kinase, sample, score) %>%
      pivot_wider(values_from = "score", names_from = "sample") %>%
      column_to_rownames("kinase")
  })

  names(act_scores) <- methods

  map(names(act_scores), function(method){
    ## Get scores for each method ---------------------------
    # change direction to perturbation
    if (method == "number_of_targets"){
      df <- map_dfr(target_df$sample, function(experiment){
        mat_meth <- data.frame(act_scores[[method]])
        colnames(mat_meth) <- colnames(act_scores[[method]])
        if (experiment %in% colnames(mat_meth)){
          mat <- t(act_scores[[method]][experiment]) * 1
          mat[mat == 0] <- NA
          data.frame(mat) %>%
            add_column(experiment = experiment, .before = 1)
        } else {
          df <- data.frame(matrix(NA, ncol = nrow(mat_meth) + 1))
          colnames(df) <- c("experiment", rownames(mat_meth))
          df
        }

      })
    } else {
      df <- map_dfr(target_df$sample, function(experiment){
        mat_meth <- data.frame(act_scores[[method]])
        colnames(mat_meth) <- colnames(act_scores[[method]])
        if (experiment %in% colnames(mat_meth)){
          mat <- t(mat_meth[experiment]) * target_df$sign[target_df$sample == experiment]

          data.frame(mat) %>%
            add_column(experiment = experiment, .before = 1)
        } else {
          df <- data.frame(matrix(NA, ncol = nrow(mat_meth) + 1))
          colnames(df) <- c("experiment", rownames(mat_meth))
          df
        }

      })
    } %>% filter(!is.na(experiment))

    df <- df[df$experiment %in% target_df$sample,]

    # filter out experiments where no activity was inferred for perturbed kinases
    if (rm_experiments){
      df_filtered <- map_dfr(1:nrow(df), function(i){
        tmp <- df[i,]
        targets <- target_df %>%
          filter(sample %in% tmp$experiment) %>%
          pull(perturb) %>%
          str_split(";") %>%
          unlist

        sum_act <- sum(tmp[,colnames(tmp) %in% targets], na.rm = T)

        if(!is.na(sum_act) & !sum_act == 0){
          df[i,]
        } else {
          df[i,] %>% mutate(experiment = "remove")
        }
      })

      df_filtered <- df_filtered %>%
        filter(!experiment == "remove")
    } else {
      df_filtered <- df
    }


    write_csv(df_filtered, paste0("results/hijazi/05_benchmark_files/subset/", method, "-", prior_i, ".csv"))

    # filter out meta to fit to experiments in matrix
    target_df <- target_df %>%
      filter(sample %in% df_filtered$experiment) %>%
      mutate(sign = 1)

    write_csv(target_df, paste0("results/hijazi/05_benchmark_files/subset/obs_", method, "-", prior_i, ".csv"))
  })

})
