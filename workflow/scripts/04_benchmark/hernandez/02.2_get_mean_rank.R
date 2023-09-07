#' Get mean/median rank for each method

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  prior_ov_file <- snakemake@input$overview
  output_file <- snakemake@output$output
  performance_per_exp <- snakemake@output$per_exp
  performance_per_kin <- snakemake@output$per_kin
  mean_rank_pdf <- snakemake@output$pdf
  bp_rank_pdf <- snakemake@output$boxplot
  methods_msk <- snakemake@params$mth
}else{
  input_files <- list.files("results/hernandez/final_scores/scaled", full.names = T)
  prior_ov_file <- "results/hernandez/overview_priors/coverage.csv"
  meta_file <- "results/hernandez/processed_data/benchmark_metadata.csv"
  output_file <- "results/hernandez/benchmark_mean_rank/mean_rank.csv"
  performance_per_exp <- "results/hernandez/benchmark_mean_rank/performance_per_exp.csv"
  performance_per_kin <- "results/hernandez/benchmark_mean_rank/performance_per_kin.csv"
  mean_rank_pdf <- "results/hernandez/benchmark_mean_rank/mean_rank.pdf"
  bp_rank_pdf <- "results/hernandez/benchmark_mean_rank/bp_rank.pdf"
  methods_msk <- c("INKA", "KARP", "KS", "KSEA_z", "PC1", "RoKAI_z", "UQ", "Wilcox", "fgsea", "mean", "median", "mlm", "norm_fgsea","norm_wmean","number_of_targets","rokai_lm","ulm","viper","wsum")
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

## Load  meta ---------------------------
obs <- read_csv(meta_file) %>%
  group_by(id) %>%
  summarise(Target = paste(target, collapse = ";"), sign = unique(sign))

## Load  activity scores ---------------------------
coverage_priors <- read_csv(prior_ov_file)
ranks <- map_dfr(input_files, function(input_file){
  net <- str_remove(str_split(input_file, "/")[[1]][5], ".rds")
  act_scores <- readRDS(input_file)
    map_dfr(names(act_scores), function(meth){
      method_act <- act_scores[[meth]]

      method_act_long <- method_act %>%
        as.data.frame() %>%
        rownames_to_column("kinase") %>%
        pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
        filter(!is.na(score)) %>%
        left_join(obs %>%
                    dplyr::rename("sample" = id) %>%
                    dplyr::select(sample, sign), by = "sample")

      if (meth == "number_of_targets"){
        method_act_long <- method_act_long %>%
          arrange(desc(score))
      } else {
        method_act_long <- method_act_long %>%
          dplyr::rename("score_raw" = score) %>%
          mutate(score = sign*score_raw) %>%
          arrange(desc(score))
      }


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
                     kinases_act = nrow(act_df)) %>%
            mutate(scaled_rank = rank/kinases_act)
        }

      })

    })
})

mean_rank_df <- ranks %>%
  filter(method %in% methods_msk) %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(method, prior) %>%
  summarise(mean_rank = round(mean(rank), digits = 1),
            mean_scaled_rank = round(mean(scaled_rank), digits = 2),
            sd_rank = round(sd(rank), digits = 1),
            sd_scaled_rank = round(sd(scaled_rank), digits = 2)) %>%
  arrange(mean_scaled_rank, mean_rank) %>%
  filter(!is.na(prior)) %>%
  mutate(comb = paste(method, prior, sep = ":"))

mean_rank_df$comb = factor(mean_rank_df$comb, levels = unique(mean_rank_df$comb))
mean_rank_df$method = factor(mean_rank_df$method, levels = unique(mean_rank_df$method))
mean_rank_df$prior = factor(mean_rank_df$prior, levels = rev(unique(mean_rank_df$prior)))

p1 <- ggplot(mean_rank_df, aes(x = mean_scaled_rank, y = method, color = prior)) +
  geom_point() +
  geom_errorbar(aes(xmin=mean_scaled_rank-sd_scaled_rank, xmax=mean_scaled_rank+sd_scaled_rank), width=.2,
                position=position_dodge(0.05)) +
  theme_minimal() +
  xlim(-0.1, 1) + facet_grid(prior ~ .)

all_ranks <- ranks %>%
  filter(method %in% methods_msk) %>%
  filter(!is.na(scaled_rank))
all_ranks$method = factor(all_ranks$method, levels = unique(mean_rank_df$method))
all_ranks$prior = factor(all_ranks$prior, levels = rev(unique(mean_rank_df$prior)))

p2 <- ggplot(all_ranks, aes(x = scaled_rank, y = method, color = prior)) +
  geom_boxplot() +
  theme_minimal() +
  xlim(-0.1, 1) + facet_grid(prior ~ .)

write_csv(mean_rank_df, output_file)

pdf(file=mean_rank_pdf, height = 17, width = 10)
plot(p1)
dev.off()

pdf(file=bp_rank_pdf, height = 17, width = 10)
plot(p2)
dev.off()

# Check experiments that perform consistently bad across priors + methods
mean_rank_exp <- ranks %>%
  filter(method %in% methods_msk) %>%
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
  filter(method %in% methods_msk) %>%
  filter(!method == "number_of_targets") %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(targets) %>%
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


