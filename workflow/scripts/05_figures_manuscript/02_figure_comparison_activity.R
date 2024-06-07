#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  benchmark_file <- snakemake@input$bench
  benchmark_meta <- snakemake@input$benchMeta
  hijazi_file <- snakemake@input$hijazi
  hijazi_meta <- snakemake@input$hijMeta
  priors <- snakemake@input$prior_files
  hijazi_priors <- snakemake@input$hijPrior
  input_file <- snakemake@input$act
  methods <- snakemake@params$methods
  overview_data_p <- snakemake@output$overview
  cor_methods_p <- snakemake@output$corrMeth
  cor_priors_p <- snakemake@output$corrPrior
}else{
  benchmark_file <- "results/01_processed_data/hernandez/data/benchmark_data.csv"
  benchmark_meta <- "results/01_processed_data/hernandez/data/benchmark_metadata.csv"
  hijazi_file <- "results/01_processed_data/hijazi/data/benchmark_data.csv"
  hijazi_meta <- "results/01_processed_data/hijazi/data/benchmark_metadata.csv"
  priors <- list.files("results/01_processed_data/hernandez/mapped_priors", pattern = ".tsv", full.names = T)
  hijazi_priors <- list.files("results/01_processed_data/hijazi/mapped_priors", pattern = ".tsv", full.names = T)
  input_file <- list.files("results/02_activity_scores/merged/final_scores", full.names = T)
  methods <- c("KARP", "KS", "KSEA_z", "PC1", "RoKAI_z", "Wilcox", "fgsea", "mean", "median", "mlm" ,"norm_wmean","number_of_targets","rokai_lm","ulm","viper","wsum", "ptmsea")
  overview_data_p <- "results/manuscript_figures/figure_2/overview_experiment.pdf"
  cor_methods_p <- "results/manuscript_figures/figure_2/corrplot_methods.pdf"
  cor_priors_p <- "results/manuscript_figures/figure_2/corrplot_priors.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)

## Figure benchmark ---------------------------
## Load data ---------------------------
hernandez <- read_csv(benchmark_file, col_types = cols()) %>%
  column_to_rownames("ID")

hijazi <- read_csv(hijazi_file, col_types = cols()) %>%
  column_to_rownames("ID")

hernandez_meta <- read_csv(benchmark_meta, col_types = cols())
hijazi_meta <- read_csv(hijazi_meta, col_types = cols())

all_targets <- map(priors, function(p_file){
  read_tsv(p_file, col_types = cols()) %>% pull(target)
}) %>%
  unlist %>%
  unique() %>%
  str_remove("\\|auto")

all_targets_hijazi <- map(hijazi_priors, function(p_file){
  read_tsv(p_file, col_types = cols()) %>% pull(target)
}) %>%
  unlist %>%
  unique() %>%
  str_remove("\\|auto")

hernandez_filtered <- hernandez[colnames(hernandez) %in% hernandez_meta$id]
hijazi_filtered <- hijazi[colnames(hijazi) %in% hijazi_meta$id]
hernandez_meta <- hernandez_meta %>%
  filter(id %in% colnames(hernandez_filtered))
hijazi_meta <- hijazi_meta %>%
  filter(id %in% colnames(hijazi_filtered))

### Coverage ---------------------------
coverage_df <- data.frame(experiment = c(colnames(hernandez_filtered),
                                         colnames(hijazi_filtered)),
                          measures_pps = c(colSums(!is.na(hernandez_filtered)),
                                           colSums(!is.na(hijazi_filtered))))

measured_in_prior <- map_dbl(colnames(hernandez_filtered), function(col_i){
  msk <- !is.na(hernandez_filtered[[col_i]])
  sum(rownames(hernandez_filtered)[msk] %in% all_targets)
})

measured_in_prior_hijazi <- map_dbl(colnames(hijazi_filtered), function(col_i){
  msk <- !is.na(hijazi_filtered[[col_i]])
  sum(rownames(hijazi_filtered)[msk] %in% all_targets_hijazi)
})

coverage_df <- coverage_df %>%
  add_column(pps_in_prior = c(measured_in_prior,
                              measured_in_prior_hijazi)) %>%
  mutate(not_in_prior = measures_pps - pps_in_prior)

coverage_df <- coverage_df %>%
  pivot_longer(!experiment, names_to = "pps_covered", values_to = "n_pps") %>%
  mutate(pps_covered = recode(pps_covered,
                              "measures_pps" = "total",
                              "pps_in_prior" = "known",
                              "not_in_prior" = "unknown"))

mean_measured_covered <- mean(coverage_df %>% filter(pps_covered == "known") %>% pull(n_pps))
sd_measured_covered <- sd(coverage_df %>% filter(pps_covered == "known") %>% pull(n_pps))

coverage_df$pps_covered <- factor(coverage_df$pps_covered, levels = c("unknown", "known"))

length(unique(coverage_df$experiment))
coverage_df %>% filter(is.na(pps_covered)) %>% pull(n_pps) %>% mean()
coverage_df %>%
  mutate(pps_covered = case_when(
    is.na(pps_covered) ~ "all",
    !is.na(pps_covered) ~ pps_covered
  )) %>%
  filter(!pps_covered == "unknown") %>%
  pivot_wider(names_from = "pps_covered", values_from = "n_pps") %>%
  mutate(percentage = known/all) %>%
  pull(percentage) %>%
  mean

overview_p <- ggplot(coverage_df %>%
                       filter(!pps_covered == "total"), aes(fill=pps_covered, y=n_pps, x=experiment)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#D3D3D3", "#6a6a6a")) +
  geom_hline(yintercept = mean_measured_covered) +
  xlab(paste0(length(unique(coverage_df$experiment)), " experiments")) +
  ylab("# of unique phosphorylation sites") +
  guides(fill=guide_legend(title="upstream kinase"))+
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        text = element_text(size = 10))

pdf(overview_data_p, width = 6, height = 3.0)
overview_p
dev.off()

## Prior comparison ---------------------------
methods <- methods[!methods == "number_of_targets"]

act_list <- map(input_file, readRDS)
names(act_list) <- map_chr(str_split(input_file, "/"), 5) %>% str_remove(".rds")

method_comparison <- map(names(act_list), function(prior_idx){
  scores <- act_list[[prior_idx]]
  scores <- scores[names(scores) %in% methods]

  scores_all <- map(names(scores), function(method_idx){
    score_long <- scores[[method_idx]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      mutate(id = paste(kinase, sample, sep = ":")) %>%
      dplyr::select(id, score)

    colnames(score_long)[2] <- method_idx
    score_long
  }) %>%
    purrr::reduce(full_join, by = "id") %>%
    column_to_rownames("id")

  cor(scores_all, use = "complete.obs")
})

names(method_comparison) <- names(act_list)

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- colnames(method_comparison[[1]])
rownames(mean_correlation) <- rownames(method_comparison[[1]])

pdf(cor_methods_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r')
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

tmp2 <- mean_correlation[colnames(mean_correlation) == "KARP"]
tmp2[!tmp2 == 1] %>%
  range()

## Prior comparison ---------------------------
prior_comparison <- map(names(act_list[[1]]), function(method_idx){
  scores <- map(act_list, function(x) x[[method_idx]])

  scores_all <- map(names(scores), function(method_idx){
    score_long <- scores[[method_idx]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      mutate(id = paste(kinase, sample, sep = ":")) %>%
      dplyr::select(id, score)

    colnames(score_long)[2] <- method_idx
    score_long
  }) %>%
    purrr::reduce(full_join, by = "id") %>%
    column_to_rownames("id")

  cor(scores_all, use = "complete.obs")
})

names(prior_comparison) <- names(act_list[[1]])


Y <- do.call(cbind, prior_comparison)
Y_array <- array(Y, dim=c(dim(prior_comparison[[1]]), length(prior_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- colnames(prior_comparison[[1]])
rownames(mean_correlation) <- rownames(prior_comparison[[1]])

pdf(cor_priors_p, height = 3, width = 3.5)
corrplot(mean_correlation, col.lim=c(0, 1), is.corr = FALSE, type = 'lower', tl.col = 'black', order = 'hclust',  cl.pos = 'r', col = COL1('Blues'))
dev.off()

mean_correlation[rownames(mean_correlation) %in% c("GPS", "phosphositeplus", "ptmsigdb"), colnames(mean_correlation) %in% c("GPS", "phosphositeplus", "ptmsigdb")]
mean_correlation
