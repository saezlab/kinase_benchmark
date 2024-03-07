#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  overview_file <- snakemake@output$ov
  rank_file <- snakemake@output$rank
  kin_rank_file <- snakemake@output$kin
  exp_rank_file <- snakemake@output$exp
}else{
  bench_files <- list.files("results/03_benchmark/hernandez/02_mean_rank",
                            recursive = TRUE, full.names = T)
  overview_file <- "results/03_benchmark/hernandez/03_benchmark_comp/overview.csv"
  rank_file <- "results/03_benchmark/hernandez/03_benchmark_comp/combined_ranks.csv"
  kin_rank_file <- "results/03_benchmark/hernandez/03_benchmark_comp/kinase_ranks.csv"
  exp_rank_file <- "results/03_benchmark/hernandez/03_benchmark_comp/experiment_ranks.csv"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
ranks <- map_dfr(bench_files, function(file){
  read_csv(file, col_types = cols())
})

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

write_csv(covered_kinases, overview_file)


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

write_csv(mean_rank_df, rank_file)


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

write_csv(mean_rank_exp, exp_rank_file)

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
  group_by(targets) %>%
  summarise(mean_rank = round(mean(rank), digits = 1),
            mean_scaled_rank = round(mean(scaled_rank), digits = 2),
            exp = paste(unique(sample), collapse = ";")) %>%
  arrange(desc(mean_scaled_rank))

write_csv(mean_rank_kin, kin_rank_file)

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
