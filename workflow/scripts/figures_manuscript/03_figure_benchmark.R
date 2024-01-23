#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench_files
  meta_file <- snakemake@input$meta
  meta_hijazi <- snakemake@input$hijazi
  auroc_plot <- snakemake@output$auroc
  overview_meta <- snakemake@output$meta_over
  rank_file <- snakemake@input$rank
  rank_plot <- snakemake@output$rankPlt
  kinase_rank <- snakemake@input$kin
  kin_file <- snakemake@output$rankKin
  prior_size <- snakemake@input$prior
}else{
  bench_files <- list.files("results/hijazi/06_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "merged")]
  auroc_plot <- "results/manuscript_figures/figure_3/auroc_res.pdf"
  prior_size <- "results/hernandez/overview_priors/coverage.csv"
  kinase_rank <- "results/hijazi/06_mean_rank/performance_per_kin_merged.csv"
  meta_file <- "results/hernandez/processed_data/benchmark_metadata.csv"
  meta_hijazi <- "results/hijazi/01_processed_data/benchmark_metadataPrior.csv"
  overview_meta <- "results/manuscript_figures/figure_3/overview_kin.pdf"
  rank_file <- "results/hijazi/06_mean_rank/full_rank_merged.csv"
  rank_plot <- "results/manuscript_figures/figure_3/mean_rank.pdf"
  kin_file <- "results/manuscript_figures/figure_3/kinase_GSknown.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)

## Overview kinases ---------------------------
if (any(str_detect(bench_files, "merged"))){
  meta_hi <- read_csv(meta_hijazi, col_types = cols()) %>%
    dplyr::select(id, sign, target)
  meta_hernandez <- read_csv(meta_file, col_types = cols()) %>%
    dplyr::select(id, sign, target)
  hernandez_meta <- rbind(meta_hi, meta_hernandez) %>%
    separate_rows(target, sep = ";")
} else if (any(str_detect(bench_files, "scaled"))){
  hernandez_meta <- read_csv(meta_file, col_types = cols()) %>%
    dplyr::select(id, sign, target)
}



kin_df <- table(hernandez_meta$target, hernandez_meta$sign) %>%
  as.data.frame() %>%
  rename("kinase" = Var1) %>%
  rename("perturbation" = Var2) %>%
  mutate(perturbation = recode(perturbation,
                               "1" = "up",
                               "-1" = "down"))
kin_order <- kin_df %>%
  group_by(kinase) %>%
  summarise(total = sum(Freq)) %>%
  arrange(desc(total)) %>%
  pull(kinase)
kin_df$kinase <- factor(kin_df$kinase, levels = rev(kin_order))

kin_p <- ggplot(kin_df, aes(x = kinase, y = Freq, fill = perturbation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  xlab("") +
  ylab("# perturbations") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))

pdf(overview_meta, width = 2.7, height = 7)
kin_p
dev.off()


## AUROC results ---------------------------
bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][4])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)

order_m <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(net)
bench_df$net <- factor(bench_df$net, levels = order_m)

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = order_method)

custom_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494")

auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  ylim(0.5, 0.83) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("AUROC") +
  xlab("")

pdf(auroc_plot, width = 6.5, height = 3)
auroc_p
dev.off()

## Mean rank ---------------------------
priors <- map_chr(str_split(bench_files, "/"), 4) %>% unique()
rank_df <- read_csv(rank_file, col_types = cols()) %>%
  filter(prior %in% priors)
rank_df$prior <- factor(rank_df$prior, levels = order_m)
rank_df$method <- factor(rank_df$method, levels = order_method)

rank_p <- rank_df %>%
  filter(!is.na(scaled_rank)) %>%
  ggplot(aes(x=method, y=scaled_rank, fill=prior)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("scaled rank") +
  xlab("")

pdf(rank_plot, width = 6.5, height = 3)
rank_p
dev.off()


## Mean rank per kinase ---------------------------
kin_rank <- read_csv(kinase_rank, col_types = cols())

kin_gsknown <- kin_rank %>%
  filter(prior == "GSknown") %>%
  arrange(mean_rank) %>%
  mutate(mean_rank = round(mean_rank)) %>%
  dplyr::select(targets, mean_rank) %>%
  mutate(range = case_when(
    mean_rank <= 10 ~ 6,
    mean_rank > 10 & mean_rank <= 20 ~ 5,
    mean_rank > 20 & mean_rank <= 30 ~ 4,
    mean_rank > 30 & mean_rank <= 40 ~ 3,
    mean_rank > 40 & mean_rank <= 50 ~ 2,
    mean_rank > 50 ~ 1
))

write_csv(kin_gsknown, kin_file)
