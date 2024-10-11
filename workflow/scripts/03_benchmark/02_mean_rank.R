#' Get mean/median rank for each method

## Snakemake ---------------------------
if(exists("snakemake")){
  input_files <- snakemake@input$scores
  target_files <- snakemake@input$target
  kinase_file <- snakemake@input$kinclass
  meta_input <- snakemake@input$meta
  rank_out <- snakemake@output$output
}else{
  input_files <- "results/03_benchmark/hernandez/01_input_bench/KS-GPS.csv"
  target_files <- "results/02_activity_scores/hernandez/scores/GPS.rds"
  kinase_file <- "resources/kinase_class.csv"
  meta_input <- "results/03_benchmark/hernandez/01_input_bench/obs_KS-GPS.csv"
  rank_out <- "results/03_benchmark/hernandez/02_mean_rank/GPS/KS-GPS.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Load scores and meta ---------------------------
method_act <- read_csv(input_files, col_types = cols()) %>%
  column_to_rownames("experiment")
targets <- readRDS(target_files)$number_of_targets %>%
  rownames_to_column("source") %>%
  pivot_longer(-source, names_to = "condition", values_to = "score") %>%
  dplyr::rename("targets" = source, "sample" = condition, "measured_targets" = score)

obs <- read_csv(meta_input, col_types = cols())

if (any(str_detect(input_files, "subset"))){
  id_col <- 6
} else {
  id_col <- 5
}

net <- str_split(str_remove(str_split(input_files, "/")[[1]][id_col], ".csv"), "-")[[1]][2]
meth <- str_split(str_remove(str_split(input_files, "/")[[1]][id_col], ".csv"), "-")[[1]][1]

kinase_class <- read_csv(kinase_file, col_types = cols()) %>%
  dplyr::filter(resource == net) %>%
  dplyr::select(source, class, kinase, n_targets) %>%
  dplyr::rename("targets" = source, "resource_class" = kinase)
  
## Get rank ---------------------------
method_act_long <- method_act %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "kinase", values_to = "score") %>%
  filter(!is.na(score))

rank_df <- map_dfr(unique(method_act_long$sample), function(exp){
  act_df <- method_act_long %>%
    dplyr::filter(sample == exp) %>%
    arrange(desc(score))

  if (exp %in% obs$sample){
    targets <- obs %>%
      dplyr::filter(sample == exp) %>%
      pull(perturb) %>%
      str_split(";") %>%
      unlist()


    rank <- map_dbl(targets, function(target){
      position <- which(act_df$kinase %in% target)
      if (length(position) == 0){
        position <- NA
      }
      position
    })

    data.frame(sample = exp,
               method = meth,
               prior = net,
               targets = targets,
               rank = rank,
               kinases_act = nrow(act_df),
               all_kinases_act = paste(act_df %>%
                                         pull(kinase), collapse = ";")) %>%
      mutate(scaled_rank = rank/kinases_act)
  }
})

rank_targets_df <- left_join(rank_df, targets, by = c("sample", "targets")) %>%
  mutate(target_size = case_when(
    measured_targets < 10 ~ "small",
    measured_targets >= 10 & measured_targets < 20 ~ "medium",
    measured_targets >= 20 ~ "large"
  ))

rank_targets_df <- left_join(rank_targets_df, kinase_class, by = c("targets"))

## Get rank ---------------------------
write_csv(rank_targets_df, rank_out)
