if(exists("snakemake")){
  meta_file <- snakemake@input$meta
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  ppsp_hernandez_file <- snakemake@input$ppspHer
  hernandez_file <- snakemake@input$her
  ppsp_hijazi_file <- snakemake@input$ppspHij
  hijazi_file <- snakemake@input$hij
  ppsp_tyrosine_file <- snakemake@input$ppspTyr
  tyrosine_file <- snakemake@input$tyr
  overview_meta <- snakemake@output$ove
  overview_meta_filtered <- snakemake@output$oveFil
  out_plot <- snakemake@output$out
}else{
  meta_file <- "results/01_processed_data/merged2/data/benchmark_metadata.csv"
  overview_meta <- "results/manuscript_figures/figure_1/supp/overview_kin.pdf"
  overview_meta_filtered <- "results/manuscript_figures/figure_1/supp/overview_kin_filtered.pdf"
  ppsp_hernandez_file <- "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv"
  hernandez_file <- "results/01_processed_data/hernandez/data/benchmark_data.csv"
  ppsp_hijazi_file <- "results/01_processed_data/hijazi/mapped_priors/phosphositeplus.tsv"
  hijazi_file <- "results/01_processed_data/hijazi/data/benchmark_data.csv"
  ppsp_tyrosine_file <- "results/01_processed_data/tyrosine/mapped_priors/phosphositeplus.tsv"
  tyrosine_file <- "results/01_processed_data/tyrosine/data/benchmark_data.csv"

  bench_files <- list.files("results/03_benchmark/hijazi/02_benchmark_res/phosphositeplus",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- c(bench_files, list.files("results/03_benchmark/hijaziDiscoverX/02_benchmark_res/phosphositeplus",
                            pattern = "bench", recursive = TRUE, full.names = T))
  rank_files <- list.files("results/03_benchmark/hijazi/02_mean_rank/phosphositeplus",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- c(rank_files, list.files("results/03_benchmark/hijaziDiscoverX/02_mean_rank/phosphositeplus",
                           pattern = "csv", recursive = TRUE, full.names = T))
  out_plot <- "results/manuscript_figures/figure_1/supp/hijazi.pdf"
}


## Libraries ---------------------------
library(tidyverse)
library(ggpubr)

## Overview kinases ---------------------------
hernandez_meta <- read_csv(meta_file, col_types = cols()) %>%
    dplyr::select(id, sign, target)
hernandez_meta <- separate_rows(hernandez_meta, target, sep = ";")

kin_df <- table(hernandez_meta$target, hernandez_meta$sign) %>%
  as.data.frame() %>%
  dplyr::rename("kinase" = Var1) %>%
  dplyr::rename("perturbation" = Var2) %>%
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
  scale_fill_manual(values=c("#E9B133", "#4E79A7")) +
  xlab("") +
  ylab("# perturbations") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))

length(unique(hernandez_meta$id))
length(unique(kin_df$kinase))
kin_df %>% group_by(perturbation) %>% summarise(total = sum(Freq))

pdf(overview_meta, width = 3, height = 8)
kin_p
dev.off()

## Overview kinases phosphositeplus ---------------------------
ppsp_hernandez <- read_tsv(ppsp_hernandez_file, col_types = cols())
ppsp_hijazi <- read_tsv(ppsp_hijazi_file, col_types = cols())
ppsp_tyrosine <- read_tsv(ppsp_tyrosine_file, col_types = cols())

data_hernandez <- read_csv(hernandez_file, col_types = cols()) %>%
  column_to_rownames("ID")
data_hijazi <- read_csv(hijazi_file, col_types = cols()) %>%
  column_to_rownames("ID")
data_tyrosine <- read_csv(tyrosine_file, col_types = cols()) %>%
  column_to_rownames("ID")

filtered_overview <- map_dfr(unique(hernandez_meta$id), function(exp){
  if (exp %in% colnames(data_hernandez)){
    data_bench <- data_hernandez
    ppsp <- ppsp_hernandez
  } else if (exp %in% colnames(data_hijazi)){    
    data_bench <- data_hijazi
    ppsp <- ppsp_hijazi
  } else if (exp %in% colnames(data_tyrosine)){
    data_bench <- data_tyrosine
    ppsp <- ppsp_tyrosine
  }
  targets <- hernandez_meta %>%
    filter(id == exp) %>% 
    pull(target)
  mat <- data_bench[exp] %>%
    drop_na()

  coverage <- map_dbl(targets, function(target_id){
    ppsp_exp <- ppsp %>%
      filter(target %in% rownames(mat))
    
    ppsp_exp %>%
      filter(source == target_id) %>%
      nrow()
  })

  data.frame(id = exp, target = targets, measure_pps = coverage)
})

hernandez_meta_filtered <- hernandez_meta %>%
  left_join(filtered_overview, by = c("id", "target")) %>%
  filter(measure_pps >= 5)

kin_df_filtered <- table(hernandez_meta_filtered$target, hernandez_meta_filtered$sign) %>%
  as.data.frame() %>%
  dplyr::rename("kinase" = Var1) %>%
  dplyr::rename("perturbation" = Var2) %>%
  mutate(perturbation = recode(perturbation,
                               "1" = "up",
                               "-1" = "down"))
kin_order <- kin_df_filtered %>%
  group_by(kinase) %>%
  summarise(total = sum(Freq)) %>%
  arrange(desc(total)) %>%
  pull(kinase)
kin_df_filtered$kinase <- factor(kin_df_filtered$kinase, levels = rev(kin_order))

kin_p_filtered <- ggplot(kin_df_filtered, aes(x = kinase, y = Freq, fill = perturbation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=c("#E9B133", "#4E79A7")) +
  xlab("") +
  ylab("# perturbations") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))

length(unique(hernandez_meta_filtered$id))
length(unique(kin_df_filtered$kinase))
kin_df_filtered %>% group_by(perturbation) %>% summarise(total = sum(Freq))

pdf(overview_meta_filtered, width = 3, height = 8)
kin_p_filtered
dev.off()

## Benchmark ------------------
## Load AUROC ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols()) %>%
    add_column(bench = str_split(file, "\\/")[[1]][3])
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets, bench) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior, "phosphositeplus" = "PhosphoSitePlus"))

n_kinases <- rank_df %>%
  group_by(bench, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  select(method, bench, kinases) %>%
  mutate(method = recode(method,
                         "chisq" = "X\u00B2 test",
                         "fisher" = "Fisher",
                         "KS" = "KS test",
                         "lmRoKAI" = "lm RoKAI",
                         "norm_mean" = "norm mean",
                         "number_of_targets" = "n targets",
                         "ptmsea" = "PTM-SEA",
                         "viper" = "VIPER",
                         "wilcox" = "MWU test",
                         "zscore" = "z-score"))

bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][5]) %>%
    add_column(bench = str_split(file, "\\/")[[1]][3])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
df_perturb_all <- bind_rows(bench_list)  %>%
  mutate(net = recode(net,"phosphositeplus" = "PhosphoSitePlus"))

df_perturb_all %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))
df_perturb_all %>% group_by(bench) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_perturb <- df_perturb_all %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(method, bench, score) %>%
  mutate(method = recode(method,
                         "chisq" = "X\u00B2 test",
                         "fisher" = "Fisher",
                         "KS" = "KS test",
                         "lmRoKAI" = "lm RoKAI",
                         "norm_mean" = "norm mean",
                         "number_of_targets" = "n targets",
                         "ptmsea" = "PTM-SEA",
                         "viper" = "VIPER",
                         "wilcox" = "MWU test",
                         "zscore" = "z-score"))

mean_auroc <- df_perturb %>%
  group_by(method, bench) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc))  # Sorting methods by mean AUROC

# Update df to ensure methods are ordered based on mean AUROC
df_perturb <- df_perturb %>%
  mutate(method = factor(method, levels = mean_auroc$method))

# Create the boxplot with ggplot2
lines <- mean_auroc$mean_auroc
names(lines) <- mean_auroc$method

auroc_p <- ggplot(df_perturb, aes(x = method, y = score, fill = bench)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5) +
  scale_fill_manual(values = c("#4292C6", "#AA42C6")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)


kin_df <- n_kinases
kin_df$method <- factor(kin_df$method, levels = mean_auroc$method)

kin_df %>% group_by(bench) %>% summarise(targets = mean(kinases))

kin_p <- ggplot(kin_df, aes(x = method, y = kinases, fill = bench)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 50),
    breaks = seq(0, 50, by = 20) #
  ) +
  theme_bw() +
  theme(legend.key.size = unit(0.3, "cm"),
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  scale_fill_manual(values = c("#4292C6", "#AA42C6")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.3, 9))

pdf(out_plot, width = 4, height = 3.8)
full_p
dev.off()
