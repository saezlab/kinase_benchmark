#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  output_file <- snakemake@output$output
  meta_out <- snakemake@output$meta_out
  rm_experiments <- snakemake@params$rm_exp
  scale_data <- snakemake@params$scale
}else{
  input_file <- "results/02_activity_scores/hernandez/final_scores/omnipath.rds"
  meta_file <- "results/01_processed_data/hernandez/data/benchmark_metadata.csv"
  meta_out <- "results/03_benchmark/hernandez/01_input_bench/obs_KARP-omnipath.csv"
  output_file <- "results/03_benchmark/hernandez/01_input_bench/KARP-omnipath.csv"
  rm_experiments <- "F"
  scale_data <- "T"
}

## Libraries ---------------------------
library(tidyverse)
source("workflow/scripts/03_benchmark/helper_functions.R")

## Define method and prior ---------------------------
method_tmp <- str_remove(str_split(output_file, "/")[[1]][5], ".csv")
input <- paste0("-",str_remove(str_split(input_file, "/")[[1]][5], ".rds"))
method <- gsub(input, "", method_tmp)

## Load scores and meta ---------------------------
act_scores <- readRDS(input_file)

if (any(names(act_scores) == "number_of_targets")){
  tmp <- act_scores["number_of_targets"]
  if(scale_data){
    act_scores <- map(act_scores, scale_scores)
  }
  act_scores["number_of_targets"] <- tmp
} else {
  if(scale_data){
    act_scores <- map(act_scores, scale_scores)
  }
}

# Scale data based on standard deviation
if(scale_data){
  act_scores <- map(act_scores, scale_scores)
}

# Prepare meta_data
obs <- read_csv(meta_file, col_types = cols())
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
      mat <- t(act_scores[[method]][,experiment]) * 1
      mat[mat == 0] <- NA
      df <- data.frame(mat) %>%
        add_column(experiment = experiment, .before = 1)
      colnames(df) <- c("experiment", rownames(mat_meth))
      df
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

