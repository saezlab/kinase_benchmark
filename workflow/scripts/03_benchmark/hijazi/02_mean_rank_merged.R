#' Get mean/median rank for each method

## Snakemake ---------------------------
if(exists("snakemake")){
  input_files <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  meta_hijazi <- snakemake@input$hijazi
  cov_kin_file <- snakemake@output$cov_kin
  prior_ov_file <- snakemake@input$overview
  output_file <- snakemake@output$output
  performance_per_exp <- snakemake@output$per_exp
  performance_per_kin <- snakemake@output$per_kin
  fullrank_file <- snakemake@output$full
}else{
  input_files <- list.files("results/hijazi/05_benchmark_files/merged", full.names = T, pattern = ".csv")
  prior_ov_file <- "results/hernandez/overview_priors/coverage.csv"
  meta_file <- "results/hernandez/processed_data/benchmark_metadata.csv"
  meta_hijazi <- "results/hijazi/01_processed_data/benchmark_metadataPrior.csv"
  output_file <- "results/hijazi/06_mean_rank/mean_rank_merged.csv"
  performance_per_exp <- "results/hijazi/06_mean_rank/performance_per_exp_merged.csv"
  performance_per_kin <- "results/hijazi/06_mean_rank/performance_per_kin_merged.csv"
  cov_kin_file <- "results/hijazi/06_mean_rank/overview/covered_kinases_merged.csv"
  fullrank_file <- "results/hijazi/06_mean_rank/full_rank_merged.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Load  meta ---------------------------
obs_hernandez <- read_csv(meta_file, col_types = cols()) %>%
  group_by(id) %>%
  summarise(Target = paste(target, collapse = ";"), sign = unique(sign))

obs_hijazi <- read_csv(meta_hijazi, col_types = cols()) %>%
  group_by(id) %>%
  summarise(Target = paste(target, collapse = ";"), sign = unique(sign))

obs <- rbind(obs_hernandez, obs_hijazi)

input_files <- input_files[!str_detect(input_files, "obs")]

## Load  activity scores ---------------------------
coverage_priors <- read_csv(prior_ov_file, col_types = cols())

## Get rank ---------------------------
ranks <- map_dfr(input_files, function(input_file){
  net <- str_split(str_remove(str_split(input_file, "/")[[1]][5], ".csv"), "-")[[1]][2]
  meth <- str_split(str_remove(str_split(input_file, "/")[[1]][5], ".csv"), "-")[[1]][1]
  method_act <- read_csv(input_file, col_types = cols()) %>%
    column_to_rownames("experiment")

  method_act_long <- method_act %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      pivot_longer(!sample, names_to = "kinase", values_to = "score") %>%
      filter(!is.na(score)) %>%
      left_join(obs %>%
                  dplyr::rename("sample" = id) %>%
                  dplyr::select(sample, sign), by = "sample") %>%
    arrange(desc(score))


    map_dfr(unique(method_act_long$sample), function(exp){
      act_df <- method_act_long %>%
        dplyr::filter(sample == exp)

      if (exp %in% obs$id){
        targets <- obs %>%
          dplyr::filter(id == exp) %>%
          pull(Target) %>%
          str_split(";") %>%
          unlist()


        rank <- map_dbl(targets, function(target){
          position <- which(act_df$kinase %in% target)
          if (length(position) == 0){
            position <- NA
          }
          position
        })


        kin_n <- coverage_priors %>%
          filter(PKN == net) %>%
          filter(class == "kinase") %>%
          filter(!type == "all kinases") %>%
          pull(value)

        data.frame(sample = exp,
                   method = meth,
                   prior = net,
                   targets = targets,
                   rank = rank,
                   kinases = kin_n,
                   kinases_act = nrow(act_df),
                   all_kinases_act = paste(method_act_long %>%
                                             filter(sample == exp) %>%
                                             pull(kinase), collapse = ";")) %>%
          mutate(scaled_rank = rank/kinases_act)
      }
    })
})

write_csv(ranks, fullrank_file)

mean_rank_df <- ranks %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(method, prior) %>%
  summarise(mean_rank = round(mean(rank), digits = 1),
            mean_scaled_rank = round(mean(scaled_rank), digits = 2),
            sd_rank = round(sd(rank), digits = 1),
            sd_scaled_rank = round(sd(scaled_rank), digits = 2),
            total_kinases = mean(kinases_act)) %>%
  arrange(mean_scaled_rank, mean_rank) %>%
  filter(!is.na(prior)) %>%
  mutate(comb = paste(method, prior, sep = ":"))

write_csv(mean_rank_df, output_file)

overview_kinases <- ranks %>%
  mutate(true_target = case_when(
    !is.na(rank) ~ targets,
    is.na(rank) ~ NA
  )) %>%
  group_by(method, prior, sample) %>%
  summarise(n_kinases = sum(!is.na(rank)),
            total_kinases = mean(kinases_act),
            target = paste(true_target %>% na.omit(), collapse = ";"),
            all_kinases = paste(unique(unlist(str_split(all_kinases_act, ";"))), collapse = ";"))

covered_kinases <- overview_kinases %>%
  group_by(prior, method) %>%
  mutate(target = case_when(
    !target == "" ~ target,
    target == "" ~ NA
  )) %>%
  summarise(TP = sum(n_kinases),
            n_targets = length(unique(unlist(str_split(target, ";"))) %>% na.omit()),
            targets =  paste(unique(unlist(str_split(target, ";"))) %>% na.omit(), collapse = ";"),
            TN = length(unlist(str_split(all_kinases, ";"))),
            n_kinases = length(unique(unlist(str_split(all_kinases, ";")))),
            all_kinases = paste(unique(unlist(str_split(all_kinases, ";"))), collapse = ";"))

write_csv(covered_kinases, cov_kin_file)

# Check experiments that perform consistently bad across priors + methods
mean_rank_exp <- ranks %>%
  filter(!method == "number_of_targets") %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(sample) %>%
  summarise(mean_rank = round(mean(rank), digits = 1),
            mean_scaled_rank = round(mean(scaled_rank), digits = 2),
            sd_rank = round(sd(rank), digits = 1),
            sd_scaled_rank = round(sd(scaled_rank), digits = 2),
            target = paste(unique(targets), collapse = ";")) %>%
  arrange(desc(mean_scaled_rank))

write_csv(mean_rank_exp, performance_per_exp)

mean_rank_exp %>%
  filter(mean_scaled_rank > 0.7) %>%
  pull(sample)

mean_rank_exp %>%
  filter(mean_scaled_rank > 0.5) %>%
  pull(sample)

# across kinases
mean_rank_kin <- ranks %>%
  filter(!method == "number_of_targets") %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(prior, targets) %>%
  summarise(mean_rank = round(mean(rank), digits = 1),
            mean_scaled_rank = round(mean(scaled_rank), digits = 2),
            exp = paste(unique(sample), collapse = ";")) %>%
  arrange(desc(mean_scaled_rank))

write_csv(mean_rank_kin, performance_per_kin)

mean_rank_kin %>%
  filter(mean_scaled_rank > 0.5) %>%
  pull(targets)

ranks %>%
  filter(targets == "CAMK2A") %>%
  pull(sample) %>%
  unique()

ranks %>%
  filter(targets == "RAF1") %>%
  pull(sample) %>%
  unique()


