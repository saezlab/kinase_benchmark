#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  output_file <- snakemake@output$output
  meta_out <- snakemake@output$meta_out
  rm_experiments <- snakemake@params$rm_exp
}else{
  input_file <- "results/hernandez/final_scores/scaled/phosphositeplus.rds"
  meta_file <- "results/hernandez/processed_data/benchmark_metadata.csv"
  meta_out <- "results/hernandez/benchmark_files/obs_KSEA_z-phosphositeplus.csv"
  output_file <- "results/hernandez/benchmark_files/KSEA_z-phosphositeplus.csv"
  rm_experiments <- "T"
}

## Libraries ---------------------------
library(tidyverse)

method_tmp <- str_remove(str_split(output_file, "/")[[1]][4], ".csv")
input <- paste0("-",str_remove(str_split(input_file, "/")[[1]][5], ".rds"))
method <- gsub(input, "", method_tmp)

rm_experiments <- as.logical(rm_experiments)
## Load scores and meta ---------------------------
act_scores <- readRDS(input_file)

obs <- read_csv(meta_file)

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
    mat <- t(act_scores[[method]][experiment]) * 1
    data.frame(mat) %>%
      add_column(experiment = experiment, .before = 1)
  })
} else {
  df <- map_dfr(target_df$sample, function(experiment){
    mat_meth <- data.frame(act_scores[[method]])
    colnames(mat_meth) <- colnames(act_scores[[method]])
    mat <- t(mat_meth[experiment]) * target_df$sign[target_df$sample == experiment]
    data.frame(mat) %>%
      add_column(experiment = experiment, .before = 1)
  })
}

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

