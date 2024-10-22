if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  kin_pdf <- snakemake@output$kin
  kinase_pdf <- snakemake@output$kin_heat
  edg_pdf <- snakemake@output$edges
  edge_pdf <- snakemake@output$edges_heat
}else{
  prior_files <- list.files("results/00_prior", pattern = "tsv", full.names = T)
  prior_files <- prior_files[c(1:4,6,7,9,10)]
  kin_pdf <- "results/manuscript_figures/figure_3/coverage_kin.pdf"
  kinase_pdf <- "results/manuscript_figures/figure_3/kinase_overview.pdf"
  edg_pdf <- "results/manuscript_figures/figure_3/coverage_edge.pdf"
  edge_pdf <- "results/manuscript_figures/figure_3/edge_overview.pdf"

  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[!str_detect(bench_files, "/shuffled/")]
  rank_files <- list.files("results/03_benchmark/merged/02_mean_rank",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[!str_detect(rank_files, "/shuffled/")]
  tumor_files <- "data/tumor_benchmark/activity_scores/roc_data.rds"
  performance_plot <- "results/04_exploration/all/performance_subset.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(ComplexHeatmap)
library(patchwork)

## Compare coverage ------------------
prior <- map(prior_files, function(file){read_tsv(file, col_types = cols()) %>%
    dplyr::filter(mor == 1) %>%
    filter(!str_detect(source, "-family")) %>%
    filter(!str_detect(source, "-subfamily"))})
names(prior) <- str_remove(str_remove(prior_files, "results/00_prior/"), ".tsv")

coverage <- map_dfr(names(prior), function(PKN_idx){
  PKN <- prior[[PKN_idx]]
  kinases <- table(PKN$source) %>% as.data.frame()

  data.frame(PKN = PKN_idx,
             value = c(nrow(kinases),
                       nrow(kinases %>% dplyr::filter(Freq >= 5)),
                       length(unique(PKN$target)),
                       nrow(PKN)),
             type = c("all kinases",
                      "kinases with \nat least 5 targets",
                      "none",
                      "none"),
             class = c("kinase", "kinase", "pps", "edges"))
})

coverage <- coverage %>% arrange(desc(value)) %>%
  mutate(PKN = recode(PKN,
                         "iKiPdb" = "iKiP-DB",
                         "GPS" = "GPS gold",
                         "omnipath" = "OmniPath",
                         "networkin" = "NetworKIN",
                         "phosphositeplus" = "PhosphoSitePlus",
                         "ptmsigdb" = "PTMsigDB",
                         "combined" = "Curated + OmniPath",
                         "GSknown" = "Curated"))
coverage %>%
  filter(type == "all kinases" & class == "kinase")

coverage %>%
  filter(class == "edges")

PKN_order <- coverage %>%
  dplyr::filter(class == "kinase") %>%
  dplyr::select(PKN, type, value) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  arrange(desc(`all kinases`), desc(`kinases with \nat least 5 targets`)) %>%
  pull(PKN)
coverage$PKN <- factor(coverage$PKN, levels = PKN_order)

text_size <- 10

kin_p <- ggplot(data=coverage %>% filter(class == "kinase"), aes(x=value, y=PKN, fill=type)) +
  geom_bar(stat="identity",color="black", position=position_dodge(), width=0.7) +
  scale_fill_manual(labels = c("all", "> 4 targets"), values=c('#477AA3','#97CAF3'))+
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, family = "Helvetica", size = 10),
        text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
        axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
        legend.position = "bottom",
        legend.key.size = unit(0.2, 'cm')) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal()

edges_p <- ggplot(coverage %>% filter(class == "edges")) +
  aes(x = value, y = PKN) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill = "#AD477A", width=0.7)+
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
        axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
        axis.title.x = element_text(family = "Helvetica", size = 10)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal()


## Kinase overview ------------------
resource_df <- map_dfr(names(prior), function(x){
  df <- prior[[x]]
  df %>%
    add_column(resource = x) %>%
    mutate(edge = paste(source, target, sep = "_"))
}) %>%
  mutate(resource = recode(resource,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "Curated + OmniPath",
                      "GSknown" = "Curated"))

kinase_m <- resource_df %>%
  dplyr::select(source,resource) %>%
  add_column(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = source, values_from = present, values_fill = 0) %>%
  column_to_rownames("resource")

kinase_tmp <- kinase_m[!str_detect(rownames(kinase_m), "Curated"), ]
sum(colSums(kinase_tmp) > 1)/ncol(kinase_tmp)
kinase_tmp[colSums(kinase_tmp) == 1] %>%
  rowSums()
order_kin <- colSums(kinase_tmp) %>%
  order()
kinase_m <- kinase_m[c(PKN_order), rev(order_kin)]

coverage_df <- pheatmap(kinase_m, cluster_rows = F,
                        cluster_cols = F, show_colnames = F,
                        treeheight_row = 0, color = c("white", "#477AA3"))


pdf(file=kin_pdf, height = 1.91, width = 3.5)
kin_p
dev.off()

pdf(file=edg_pdf, height = 1.91, width = 2.3)
edges_p
dev.off()

pdf(kinase_pdf, height = 1.5, width = 3.7)
print(coverage_df)
dev.off()


## Edge overview ------------------
edge_m <- resource_df %>%
  dplyr::select(edge,resource) %>%
  add_column(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = edge, values_from = present, values_fill = 0) %>%
  column_to_rownames("resource")

order_edge <- colSums(edge_m) %>%
  order()
edge_m <- edge_m[c(PKN_order), rev(order_edge)]
sum(colSums(edge_m) > 1)/ncol(edge_m)
edge_m[colSums(edge_m) == 1] %>%
  rowSums()

coverage_edge_df <- pheatmap(edge_m, cluster_rows = F,
                             cluster_cols = F, show_colnames = F,
                             treeheight_row = 0, color = c("white", "#AD477A"))

pdf(edge_pdf, height = 1.5, width = 10)
print(coverage_edge_df)
dev.off()

## Benchmark ------------------
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
                      "combined" = "Curated + OmniPath",
                      "GSknown" = "Curated",
                      "shuffled2" = "Shuffled"))


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
                        "iKiPdb" = "iKiP-DB",
                        "GPS" = "GPS gold",
                        "omnipath" = "OmniPath",
                        "networkin" = "NetworKIN",
                        "phosphositeplus" = "PhosphoSitePlus",
                        "ptmsigdb" = "PTMsigDB",
                        "combined" = "Curated + OmniPath",
                        "GSknown" = "Curated",
                        "shuffled2" = "Shuffled"))

df_perturb_all %>% filter(!net %in% c("shuffled", "shuffled2")) %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_perturb <- df_perturb_all %>%
  dplyr::filter(method == "zscore") %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(net, score) %>%
  add_column(benchmark = "perturbation-based")

known_roc <- readRDS(tumor_files)
roc_list <- map(known_roc, function(known_roc_i) {
  lapply(known_roc_i$ROC_results, "[[", "sample_AUROCs")}
)
names(roc_list) <- names(known_roc)

df_tumor_all <- map_dfr(names(roc_list), function(roc_i){
  roc <- roc_list[[roc_i]]
  do.call(rbind, lapply(names(roc), function(method) {
    data.frame(method = method, score = roc[[method]], benchmark = "tumor-based")
  })) %>%
    add_column(net = roc_i)
})
df_tumor_all %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_tumor <- df_tumor_all %>%
  dplyr::filter(method == "zscore") %>%
  dplyr::select(net, score, benchmark) %>%
  mutate(net = recode(net,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "Curated + OmniPath",
                      "GSknown" = "Curated",
                      "shuffled2" = "Shuffled"))

n_kinases_tumor <- map_dfr(names(roc_list), function(roc_i){
  roc <- known_roc[[roc_i]]
  data.frame(prior = roc_i,
             kinases = length(unique(roc$evaluation_kinases %>% unlist())),
             benchmark = "tumor-based")
}) %>%
  mutate(prior = recode(prior,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "Curated + OmniPath",
                      "GSknown" = "Curated",
                      "shuffled2" = "Shuffled"))

bench_df <- rbind(df_perturb, df_tumor)

mean_auroc <- bench_df %>%
  group_by(net, benchmark) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(net) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc))  # Sorting methods by mean AUROC

# Update df to ensure methods are ordered based on mean AUROC
bench_df <- bench_df %>%
  mutate(net = factor(net, levels = mean_auroc$net))

# Create the boxplot with ggplot2
lines <- mean_auroc$mean_auroc
names(lines) <- mean_auroc$net

auroc_p <- ggplot(bench_df, aes(x = net, y = score, fill = benchmark)) +
  geom_boxplot(linewidth = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c("#4292C6", "#b54d4a")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
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


kin_df <- rbind(n_kinases, n_kinases_tumor)
kin_df$prior <- factor(kin_df$prior, levels = mean_auroc$net)

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
  scale_fill_manual(values = c("#4292C6", "#b54d4a")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf("results/manuscript_figures/figure_3/zscore.pdf", width = 6.5, height = 4)
full_p
dev.off()
