#' Get mean/median rank for each method

## Snakemake ---------------------------
if(exists("snakemake")){
  input_files <- snakemake@input$scores
  meta_input <- snakemake@input$meta
  rank_out <- snakemake@output$output
}else{
  input_files <- "results/03_benchmark/hernandez/01_input_bench/fgsea-GPS.csv"
  meta_input <- "results/03_benchmark/hernandez/01_input_bench/obs_fgsea-GPS.csv"
  rank_out <- "results/03_benchmark/hernandez/02_mean_rank/GPS/fgsea-GPS.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Load scores and meta ---------------------------
method_act <- read_csv(input_files, col_types = cols()) %>%
  column_to_rownames("experiment")

obs <- read_csv(meta_input, col_types = cols())

net <- str_split(str_remove(str_split(input_files, "/")[[1]][5], ".csv"), "-")[[1]][2]
meth <- str_split(str_remove(str_split(input_files, "/")[[1]][5], ".csv"), "-")[[1]][1]


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


## Get rank ---------------------------
write_csv(rank_df, rank_out)
