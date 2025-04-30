if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  activating_files <- snakemake@input$act
  activating_files_kin <- snakemake@input$act_kin
  tumor_files <- snakemake@input$tumor
  tumor_files_kin <- snakemake@input$tumor_kin
  kin_class_file <- snakemake@input$kinclass
  performance_plot <- snakemake@output$plt
  overview_csv <- snakemake@output$csv
}else{
  bench_files <- list.files("results/03_benchmark/merged2/02_benchmark_res_subset/subset",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "/shuffled2/|/iKiPdb/|/GPS/|/omnipath/|/networkin/|/phosphositeplus/|/ptmsigdb/|/GSknown/")]
  bench_files <- bench_files[str_detect(bench_files, "zscore")]
  rank_files <- list.files("results/03_benchmark/merged2/02_mean_rank_subset/subset",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[str_detect(rank_files, "/shuffled2/|/iKiPdb/|/GPS/|/omnipath/|/networkin/|/phosphositeplus/|/ptmsigdb/|/GSknown/")]
  rank_files <- rank_files[str_detect(rank_files, "zscore")]
  activating_files <- list.files("data/results_cptac/overall_performance/actsiteBM/same_kins",
                                 pattern = "table", full.names = T)
  activating_files_kin <- list.files("data/results_cptac/overall_performance/actsiteBM/same_kins",
                                 pattern = "kins", full.names = T)
  tumor_files <- list.files("data/results_cptac/overall_performance/protBM/same_kins",
                                 pattern = "table", full.names = T)
  tumor_files_kin <- list.files("data/results_cptac/overall_performance/protBM/same_kins",
                                 pattern = "kins", full.names = T)
  performance_plot <- "results/manuscript_figures/figure_3/zscore.pdf"
  overview_csv <- "results/manuscript_figures/supp_files/overview_benchmark.csv"
  kin_class_file <- "resources/kinase_class.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

## Benchmark ------------------
#kin_class <- read_csv(kin_class_file, col_types = cols())
#kin_ser <- kin_class %>% filter(class == "Serine/Threonine") %>% pull(source) %>% unique()
#kin_tyr <- kin_class %>% filter(class == "Tyrosine") %>% pull(source) %>% unique()
#kin_dual <- kin_class %>% filter(class == "Dual-specificity") %>% pull(source) %>% unique()
## Load AUROC ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols())
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "extended combined",
                      "GSknown" = "curated",
                      "shuffled2" = "shuffled"))


n_kinases <- rank_df %>%
  group_by(prior, method) %>%
  summarise(TP = n(), kinases = length(unique(targets))) %>%
  filter(method == "zscore") %>%
  select(prior, kinases, TP) %>%
  add_column(benchmark = "perturbation-based")

#n_kinases_class <- rank_df %>%
#  left_join(kin_class %>% dplyr::select(source, class), by = c("targets" = "source"), relationship = "many-to-many") %>%
#  group_by(prior, method, class) %>%
#  summarise(TP = n(), kinases = length(unique(targets))) %>%
#  filter(method == "zscore") %>%
#  select(prior, kinases, TP, class) %>%
#  add_column(benchmark = "perturbation-based")

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
                        "iKiPdb" = "iKiP-DB",
                        "GPS" = "GPS gold",
                        "omnipath" = "OmniPath",
                        "networkin" = "NetworKIN",
                        "phosphositeplus" = "PhosphoSitePlus",
                        "ptmsigdb" = "PTMsigDB",
                        "combined" = "extended combined",
                        "GSknown" = "curated",
                        "shuffled2" = "shuffled"))

df_perturb_all %>% filter(!net %in% c("shuffled", "shuffled2")) %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_perturb <- df_perturb_all %>%
  dplyr::filter(method == "zscore") %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(net, score) %>%
  add_column(benchmark = "perturbation-based")

## tumor benchmark
df_tumor <- map_dfr(tumor_files, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
    data.frame(net = net_id, score = roc[,colnames(roc) == "zscore"], benchmark = "tumor-based")
}) %>%
  mutate(net = recode(net,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled"))

## load number of kinases
n_kinases_tumor <- map_dfr(tumor_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist()
  net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
  data.frame(prior = net_id,
             kinases = length(unique(kin)),
             TP = length(kin),
             benchmark = "tumor-based")
}) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled"))

#n_kinases_tumor_class <- map_dfr(tumor_files_kin, function(file_roc){
#  kin <- readRDS(file_roc) %>%
#    unlist()
#  net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
#  data.frame(prior = net_id,
#             kinases = c(sum(unique(kin) %in% kin_ser), sum(unique(kin) %in% kin_tyr), sum(unique(kin) %in% kin_dual)),
#             TP = c(sum((kin) %in% kin_ser), sum((kin) %in% kin_tyr), sum((kin) %in% kin_dual)),
#             class = c("Serine/Threonine", "Tyrosine", "Dual-specificity"),
#             benchmark = "tumor-based")
#}) %>%
#  mutate(prior = recode(prior,
#                      "ikip" = "iKiP-DB",
#                      "gps" = "GPS gold",
#                      "omni" = "OmniPath",
#                      "nwkin" = "NetworKIN",
#                      "psp" = "PhosphoSitePlus",
#                      "ptmsig" = "PTMsigDB",
#                      "combo" = "extended combined",
#                      "known" = "curated",
#                      "shuffled2" = "shuffled"))

## activating sites benchmark
df_act <- map_dfr(activating_files, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
    data.frame(net = net_id, score = roc[,colnames(roc) == "zscore"], benchmark = "activating sites")
}) %>%
  mutate(net = recode(net,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled"))

## load number of kinases
n_kinases_act <- map_dfr(activating_files_kin, function(file_roc){
  kin <- readRDS(file_roc) %>%
    unlist()
  net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
  data.frame(prior = net_id,
             kinases = length(unique(kin)),
             TP = length(kin),
             benchmark = "activating sites")
}) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled"))

#n_kinases_act_class <- map_dfr(activating_files_kin, function(file_roc){
#  kin <- readRDS(file_roc) %>%
#    unlist()
#  net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
#  data.frame(prior = net_id,
#             kinases = c(sum(unique(kin) %in% kin_ser), sum(unique(kin) %in% kin_tyr), sum(unique(kin) %in% kin_dual)),
#             TP = c(sum((kin) %in% kin_ser), sum((kin) %in% kin_tyr), sum((kin) %in% kin_dual)),
#             class = c("Serine/Threonine", "Tyrosine", "Dual-specificity"),
#             benchmark = "activating sites")
#}) %>%
#  mutate(prior = recode(prior,
#                      "ikip" = "iKiP-DB",
#                      "gps" = "GPS gold",
#                      "omni" = "OmniPath",
#                      "nwkin" = "NetworKIN",
#                      "psp" = "PhosphoSitePlus",
#                      "ptmsig" = "PTMsigDB",
#                      "combo" = "extended combined",
#                      "known" = "curated",
#                      "shuffled2" = "shuffled"))



## Combine
df_tumor <- df_tumor %>% filter(net %in% unique(df_perturb$net))
bench_df <- rbind(df_perturb, df_tumor, df_act)

bench_df %>%
  group_by(net, benchmark) %>%
  summarise(score = mean(score)) %>%
  pivot_wider(names_from = benchmark, values_from = score) %>%
  column_to_rownames("net") %>%
  cor()

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
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)

n_kinases_tumor <- n_kinases_tumor %>% filter(prior %in% unique(n_kinases$prior))
kin_df <- rbind(n_kinases, n_kinases_tumor, n_kinases_act)
kin_df$prior <- factor(kin_df$prior, levels = mean_auroc$net)
kin_df$benchmark <- factor(kin_df$benchmark, levels = c("perturbation-based", "activating sites", "tumor-based"))
print("test 4")
if (any(str_detect(activating_files_kin, "same_kins"))){
  limits <- c(0,40)
  breaks <- seq(0, 50, by = 20)
} else {
  limits <- c(0, 160)
  breaks <- seq(0, 160, by = 50) 
}

kin_p <- ggplot(kin_df, aes(x = prior, y = kinases, fill = benchmark)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = limits,
    breaks = breaks#
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


if (!str_detect(performance_plot, "supp")){
  write_csv(bench_df, "results/manuscript_data/fig3c_boxplot.csv")
  write_csv(kin_df %>% dplyr::select(-TP), "results/manuscript_data/fig3c_barplot.csv")
} else if (str_detect(performance_plot, "supp")){
  write_csv(bench_df, "results/manuscript_data/suppfig9_boxplot.csv")
  write_csv(kin_df %>% dplyr::select(-TP), "results/manuscript_data/suppfig9_barplot.csv")
} 



full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf(performance_plot, width = 3.9, height = 4.2)
full_p
dev.off()

bench_df %>% group_by(net) %>% summarise(score = mean(score)) %>% arrange(desc(score))
kin_df %>% filter(prior %in% c("curated"))
kin_df %>% filter(prior %in% c("OmniPath", "iKiP-DB"))

#kin_df_class <- rbind(n_kinases_class, n_kinases_tumor_class, n_kinases_act_class)

#if (!str_detect(performance_plot, "/supp/")){
#  kin_csv <- kin_df_class %>%
#    dplyr::select(-method) %>% 
#    dplyr::mutate(benchmark = recode(benchmark,
#                              "tumor-based" = "protein-based",
#                              "activating sites" = "activating site-based")) %>%
#    dplyr::rename("kinase-substrate library" = prior,
#          "unique kinases in GS set" = kinases,
#          "unique pairs in GS set" = TP,
#          "kinase class" = class,
#          "benchmark approach" = benchmark)
#
#  write_csv(kin_csv, overview_csv)
#}

