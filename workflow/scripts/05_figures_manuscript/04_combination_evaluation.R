if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  activating_files <- snakemake@input$act
  tumor_files <- snakemake@input$tumor
  performance_plot <- snakemake@output$plt
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "/GSknown/|/GSknown_iKiPdb/|/GSknown_networkin/|/combined/")]
  rank_files <- list.files("results/03_benchmark/merged/02_mean_rank",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[str_detect(rank_files, "/GSknown/|/GSknown_iKiPdb/|/GSknown_networkin/|/combined/")]
  activating_files <- list.files("data/results_cptac/performance_combinations/GSknown/all_kins",
                                 pattern = "actsiteBM", full.names = T)
  activating_files_roc <- activating_files[!str_detect(activating_files, "known_KL|known_phosf|kins.rds")]
  activating_files_kin <- activating_files[!str_detect(activating_files, "known_KL|known_phosf|table.rds")]
  tumor_files <- list.files("data/results_cptac/performance_combinations/GSknown/all_kins",
                                 pattern = "protBM", full.names = T)
  tumor_files_roc <- tumor_files[!str_detect(tumor_files, "known_KL|known_phosf|kins.rds")]
  tumor_files_kin <- tumor_files[!str_detect(tumor_files, "known_KL|known_phosf|table.rds")]
  performance_plot <- "results/manuscript_figures/figure_4/combinations_zscore.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

## Benchmark ------------------
## Load AUROC ---------------------------
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
                        "combined" = "curated + OmniPath"))


n_kinases <- rank_df %>%
  group_by(prior, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  filter(method == "zscore") %>%
  select(prior, kinases) %>%
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
                        "combined" = "curated + OmniPath"))

df_perturb_all %>% filter(!net %in% c("shuffled", "shuffled2")) %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_perturb <- df_perturb_all %>%
  dplyr::filter(method == "zscore") %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(net, score) %>%
  add_column(benchmark = "perturbation-based")

## tumor benchmark
df_tumor <- map_dfr(tumor_files_roc, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_table)")
    data.frame(net = net_id, score = roc[colnames(roc) == "zscore"], benchmark = "tumor-based")
}) %>%
   mutate(net = recode(net,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath"))


## load number of kinases
n_kinases_tumor <- map_dfr(tumor_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist() %>%
    unique()
  net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_kins)")
  data.frame(prior = net_id,
             kinases = length(kin),
             benchmark = "tumor-based")
}) %>%
   mutate(prior = recode(prior,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath"))


## activating sites benchmark
df_act <- map_dfr(activating_files_roc, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_table)")
    data.frame(net = net_id, score = roc[colnames(roc) == "zscore"], benchmark = "activating sites")
}) %>%
   mutate(net = recode(net,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath"))


## load number of kinases
n_kinases_act <- map_dfr(activating_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist() %>%
    unique()
  net_id <- str_extract(file_roc, "(?<=combofig_).*?(?=_roc_kins)")
  data.frame(prior = net_id,
             kinases = length(kin),
             benchmark = "activating sites")
}) %>%
   mutate(prior = recode(prior,
                        "known" = "curated",
                        "known_ikip" = "curated + iKiP-DB",
                        "known_nwkin" = "curated + NetworKIN",
                        "combo" = "curated + OmniPath"))


## Combine
bench_df <- rbind(df_perturb, df_tumor, df_act)

mean_auroc <- bench_df %>%
  group_by(net, benchmark) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(net) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc))  # Sorting methods by mean AUROC

# Update df to ensure methods are ordered based on mean AUROC
bench_df <- bench_df %>%
  mutate(net = factor(net, levels = mean_auroc$net)) %>%
  mutate(benchmark = factor(benchmark, levels = c("perturbation-based", "activating sites", "tumor-based")))

# Create the boxplot with ggplot2
lines <- mean_auroc$mean_auroc
names(lines) <- mean_auroc$net

auroc_p <- ggplot(bench_df, aes(x = net, y = score, fill = benchmark)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#4292C6", "#C67642", "#b54d4a")) +  # Muted scientific color palette
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
  # Add horizontal lines or dots for mean AUROC
  geom_point(data = data.frame(net = names(lines), mean_auroc = lines), aes(x = net, y = mean_auroc),
             color = "black", size = 2, shape = 3, fill = "white")  +
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)


kin_df <- rbind(n_kinases, n_kinases_tumor, n_kinases_act)
kin_df$prior <- factor(kin_df$prior, levels = mean_auroc$net)
kin_df$benchmark <- factor(kin_df$benchmark, levels = c("perturbation-based", "activating sites", "tumor-based"))

kin_p <- ggplot(kin_df, aes(x = prior, y = kinases, fill = benchmark)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 175),
    breaks = seq(0, 175, by = 50) #
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
  scale_fill_manual(values = c("#4292C6", "#C67642", "#b54d4a")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf(performance_plot, width = 3.9, height = 4.2)
full_p
dev.off()
