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
  input_file <- "results/hijazi/04_final_scores/scaled/GSknown.rds"
  input_hernandez <- "results/hernandez/final_scores/scaled/GSknown.rds"
  meta_file <- "results/hijazi/01_processed_data/benchmark_metadataPrior.csv"
  meta_hernandez <- "results/hernandez/processed_data/benchmark_metadata.csv"
  meta_out <- "results/hijazi/05_benchmark_files/merged/obs_number_of_targets-GSknown.csv"
  output_file <- "results/hijazi/05_benchmark_files/merged/number_of_targets-GSknown.csv"
  rm_experiments <- "F"
}

## Libraries ---------------------------
library(tidyverse)

method_tmp <- str_remove(str_split(output_file, "/")[[1]][5], ".csv")
input <- paste0("-",str_remove(str_split(input_file, "/")[[1]][5], ".rds"))
method <- gsub(input, "", method_tmp)

rm_experiments <- as.logical(rm_experiments)
## Load scores and meta ---------------------------
act_scores_hijazi <- readRDS(input_file)
act_scores_hernandez <- readRDS(input_hernandez)

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


write_csv(df_filtered, output_file)

# filter out meta to fit to experiments in matrix
target_df <- target_df %>%
  filter(sample %in% df_filtered$experiment) %>%
  mutate(sign = 1)

write_csv(target_df, meta_out)

