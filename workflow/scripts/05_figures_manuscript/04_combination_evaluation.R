if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  activating_files_roc <- snakemake@input$act_roc
  activating_files_kin <- snakemake@input$act_kin
  tumor_files_roc <- snakemake@input$tumor_roc
  tumor_files_kin <- snakemake@input$tumor_kin
  performance_plot <- snakemake@output$plt
  kin_class_file <- snakemake@input$kinclass
  performance_plot_prot <- snakemake@output$plt_prot
  performance_plot_act <- snakemake@output$plt_act
}else{
  bench_files <- list.files("results/03_benchmark/merged2/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "/GSknown/|/GSknown_johnson15/|/GSknown_phosformer15/")]
  rank_files <- list.files("results/03_benchmark/merged2/02_mean_rank",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[str_detect(rank_files, "/GSknown/|/GSknown_johnson15/|/GSknown_phosformer15/")]
  activating_files <- list.files("data/results_cptac/performance_combinations/GSknown/all_kins",
                                 pattern = "actsiteBM", full.names = T)
  activating_files_roc <- activating_files[!str_detect(activating_files, "known_ikip|known_nwkin|kins.rds")]
  activating_files_kin <- activating_files[!str_detect(activating_files, "known_ikip|known_nwkin|table.rds")]
  tumor_files <- list.files("data/results_cptac/performance_combinations/GSknown/all_kins",
                                 pattern = "protBM", full.names = T)
  tumor_files_roc <- tumor_files[!str_detect(tumor_files, "known_ikip|known_nwkin|kins.rds")]
  tumor_files_kin <- tumor_files[!str_detect(tumor_files, "known_ikip|known_nwkin|table.rds")]
  performance_plot <- "results/manuscript_figures/figure_4/combinations_zscore_perturbation.pdf"
  performance_plot_prot <- "results/manuscript_figures/figure_4/combinations_zscore_protein.pdf"
  performance_plot_act <- "results/manuscript_figures/figure_4/combinations_zscore_actsite.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

## Benchmark ------------------
kin_class <- read_csv(kin_class_file, col_types = cols())
kin_ser <- kin_class %>% filter(class == "Serine/Threonine") %>% pull(source) %>% unique()
kin_tyr <- kin_class %>% filter(class == "Tyrosine") %>% pull(source) %>% unique()
kin_dual <- kin_class %>% filter(class == "Dual-specificity") %>% pull(source) %>% unique()
## Load AUROC ---------------------------
method_selection <- c("z-score", "PTM-SEA", "mean", "VIPER", "fgsea", "KSEA", "n targets")
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols())
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior,
                        "GSknown" = "curated",
                        "GSknown_iKiPdb" = "curated + iKiP-DB",
                        "GSknown_networkin" = "curated + NetworKIN",
                        "combined" = "curated + OmniPath",
                        "GSknown_johnson15" = "curated + KinaseLibrary",
                        "GSknown_phosformer15" = "curated + Phosformer")) %>%
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



n_kinases <- rank_df %>%
  group_by(prior, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  filter(method %in% method_selection) %>%
  select(prior, method, kinases) %>%
  add_column(benchmark = "perturbation-based")

bind_rows(rank_list) %>%
  filter(!is.na(rank)) %>%
  dplyr::select(sample, prior, all_kinases_act) %>%
  separate_rows(all_kinases_act, sep = ";") %>%
  distinct() %>%
  left_join(kin_class %>% dplyr::select(source, class), by = c("all_kinases_act" = "source"), relationship = "many-to-many") %>%
  distinct() %>%
  group_by(sample, prior, class) %>%
  summarise(kinases = n()) %>%
  group_by(prior, class) %>%
  summarise(avg_kinases = mean(kinases))

bind_rows(rank_list) %>%
  filter(!is.na(rank)) %>%
  dplyr::select(prior, all_kinases_act) %>%
  separate_rows(all_kinases_act, sep = ";") %>%
  distinct() %>%
  left_join(kin_class %>% dplyr::select(source, class), by = c("all_kinases_act" = "source"), relationship = "many-to-many") %>%
  group_by(prior, class) %>%
  distinct() %>%
  summarise(kinases = n())

rank_df %>%
  left_join(kin_class %>% dplyr::select(source, class), by = c("targets" = "source"), relationship = "many-to-many") %>%
  group_by(prior, method, class) %>%
  summarise(TP = n(), kinases = length(unique(targets))) %>%
  filter(method == "z-score") %>%
  select(prior, kinases, TP, class) %>%
  add_column(benchmark = "perturbation-based")

if (any(str_detect(bench_files, "subset"))){
  net_id <- 6
} else {
  net_id <- 5
}

bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][net_id])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
df_perturb_all <- bind_rows(bench_list)  %>%
  mutate(net = recode(net,
                        "GSknown" = "curated",
                        "GSknown_iKiPdb" = "curated + iKiP-DB",
                        "GSknown_networkin" = "curated + NetworKIN",
                        "combined" = "curated + OmniPath",
                        "GSknown_johnson15" = "curated + KinaseLibrary",
                        "GSknown_phosformer15" = "curated + Phosformer")) %>%
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

df_perturb_all %>% filter(!net %in% c("shuffled", "shuffled2")) %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_perturb <- df_perturb_all %>%
  filter(method %in% method_selection) %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(net, method, score) %>%
  add_column(benchmark = "perturbation-based")



## tumor benchmark
df_tumor <- map_dfr(tumor_files_roc, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_table)")
    map_dfr(colnames(roc), function(method_id){
      data.frame(net = net_id, method = method_id, score = roc[,colnames(roc) == method_id], benchmark = "tumor-based")
    })
}) %>%
   mutate(net = recode(net,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath",
                        "known_KL" = "curated + KinaseLibrary",
                        "known_phosf" = "curated + Phosformer"))  %>%
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
                         "zscore" = "z-score")) %>%
  filter(method %in% method_selection)

## load number of kinases
n_kinases_tumor <- map_dfr(tumor_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist() %>%
    unique()
  net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_kins)")
  data.frame(prior = net_id,
             kinases = length(kin),
             benchmark = "tumor-based",
             method = unique(df_tumor$method))
}) %>%
   mutate(prior = recode(prior,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath",
                        "known_KL" = "curated + KinaseLibrary",
                        "known_phosf" = "curated + Phosformer"))

map_dfr(tumor_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist()
  net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
  data.frame(prior = net_id,
             kinases = c(sum(unique(kin) %in% kin_ser), sum(unique(kin) %in% kin_tyr), sum(unique(kin) %in% kin_dual)),
             TP = c(sum((kin) %in% kin_ser), sum((kin) %in% kin_tyr), sum((kin) %in% kin_dual)),
             class = c("Serine/Threonine", "Tyrosine", "Dual-specificity"),
             benchmark = "tumor-based")
})

## activating sites benchmark
df_act <- map_dfr(activating_files_roc, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_table)")
        map_dfr(colnames(roc), function(method_id){
      data.frame(net = net_id, method = method_id, score = roc[,colnames(roc) == method_id], benchmark = "activating sites")
    })
}) %>%
   mutate(net = recode(net,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath",
                        "known_KL" = "curated + KinaseLibrary",
                        "known_phosf" = "curated + Phosformer")) %>%
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
                         "zscore" = "z-score")) %>%
  filter(method %in% method_selection)


## load number of kinases
n_kinases_act <- map_dfr(activating_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist() %>%
    unique()
  net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_kins)")
  data.frame(prior = net_id,
             kinases = length(kin),
             benchmark = "activating sites",
             method = unique(df_act$method))
}) %>%
   mutate(prior = recode(prior,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath",
                        "known_KL" = "curated + KinaseLibrary",
                        "known_phosf" = "curated + Phosformer"))

map_dfr(activating_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist()
  net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
  data.frame(prior = net_id,
             kinases = c(sum(unique(kin) %in% kin_ser), sum(unique(kin) %in% kin_tyr), sum(unique(kin) %in% kin_dual)),
             TP = c(sum((kin) %in% kin_ser), sum((kin) %in% kin_tyr), sum((kin) %in% kin_dual)),
             class = c("Serine/Threonine", "Tyrosine", "Dual-specificity"),
             benchmark = "tumor-based")
})

## Combine
bench_df <- rbind(df_perturb, df_tumor, df_act)

mean_auroc <- bench_df %>%
  group_by(net, benchmark) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(net) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc))   %>% # Sorting methods by mean AUROC
  arrange(desc(net == "curated"))

mean_method <- bench_df %>%
  group_by(method, benchmark) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc)) 
# Update df to ensure methods are ordered based on mean AUROC
bench_df <- bench_df %>%
  mutate(net = factor(net, levels = mean_auroc$net)) %>%
  mutate(method = factor(method, levels = mean_method$method)) %>%
  mutate(benchmark = factor(benchmark, levels = c("perturbation-based", "activating sites", "tumor-based")))

# Create the boxplot with ggplot2
lines <- mean_auroc$mean_auroc
names(lines) <- mean_auroc$net

auroc_p <- ggplot(bench_df %>% filter(benchmark == "perturbation-based"), aes(x = method, y = score, fill = net)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values =  c("#ECF4F9", "#B3D3E8", "#7BB3D7", "#4292C6")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)


kin_df <- rbind(n_kinases, n_kinases_tumor, n_kinases_act)
kin_df$prior <- factor(kin_df$prior, levels = mean_auroc$net)
kin_df$method <- factor(kin_df$method, levels = mean_method$method)
kin_df$benchmark <- factor(kin_df$benchmark, levels = c("perturbation-based", "activating sites", "tumor-based"))

if("curated + Phosformer" %in% kin_df$prior){
    limits <- c(0, 80)
    breaks <- seq(0, 80, by = 20)
} else {
    limits <- c(0, 60)
    breaks <- seq(0, 60, by = 20) 
}

kin_p <- ggplot(kin_df %>% filter(benchmark == "perturbation-based"), aes(x = method, y = kinases, fill = prior)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = limits,
    breaks = breaks #
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  scale_fill_manual(values =  c("#ECF4F9", "#B3D3E8", "#7BB3D7", "#4292C6"))  +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf(performance_plot, width = 3.9, height = 3.7)
full_p
dev.off()

## activating site
auroc_p <- ggplot(bench_df %>% filter(benchmark == "tumor-based"), aes(x = method, y = score, fill = net)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values =  c("#FFECEC", "#E1B8B7", "#CB8381", "#B54D4A")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)



kin_df <- rbind(n_kinases, n_kinases_tumor, n_kinases_act)
kin_df$prior <- factor(kin_df$prior, levels = mean_auroc$net)
kin_df$method <- factor(kin_df$method, levels = mean_method$method)
kin_df$benchmark <- factor(kin_df$benchmark, levels = c("perturbation-based", "activating sites", "tumor-based"))

if("curated + Phosformer" %in% kin_df$prior){
    limits <- c(0, 300)
    breaks <- seq(0, 300, by = 100)
} else {
    limits <- c(0, 175)
    breaks <- seq(0, 175, by = 50) 
}

kin_p <- ggplot(kin_df %>% filter(benchmark == "tumor-based"), aes(x = method, y = kinases, fill = prior)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = limits,
    breaks = breaks #
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  scale_fill_manual(values =  c("#FFECEC", "#E1B8B7", "#CB8381", "#B54D4A"))   +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf(performance_plot_prot, width = 3.9, height = 3.7)
full_p
dev.off()

## protein-based
auroc_p <- ggplot(bench_df %>% filter(benchmark == "activating sites"), aes(x = method, y = score, fill = net)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values =  c("#FEE6D8", "#E8C8B3", "#D79F7B", "#C67642")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)


kin_df <- rbind(n_kinases, n_kinases_tumor, n_kinases_act)
kin_df$prior <- factor(kin_df$prior, levels = mean_auroc$net)
kin_df$method <- factor(kin_df$method, levels = mean_method$method)
kin_df$benchmark <- factor(kin_df$benchmark, levels = c("perturbation-based", "activating sites", "tumor-based"))

if("curated + Phosformer" %in% kin_df$prior){
    limits <- c(0, 90)
    breaks <- seq(0, 90, by = 30)
} else {
    limits <- c(0, 80)
    breaks <- seq(0, 80, by = 25) 
}

kin_p <- ggplot(kin_df %>% filter(benchmark == "activating sites"), aes(x = method, y = kinases, fill = prior)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = limits,
    breaks = breaks #
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  scale_fill_manual(values =  c("#FEE6D8", "#E8C8B3", "#D79F7B", "#C67642"))   +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf(performance_plot_act, width = 3.9, height = 3.7)
full_p
dev.off()
