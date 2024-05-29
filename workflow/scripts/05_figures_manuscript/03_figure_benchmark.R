#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench_files
  meta_file <- snakemake@input$meta
  auroc_plot <- snakemake@output$auroc
  overview_meta <- snakemake@output$meta_over
  rank_files <- snakemake@input$rank
  rank_plot <- snakemake@output$rankPlt
  kin_file <- snakemake@output$rankKin
  heat_plot <- snakemake@output$heat
  medRank_plot <- snakemake@output$medRank
  gene_citations <- snakemake@input$cit
  cor_plot <- snakemake@output$cor
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  meta_file <- "results/01_processed_data/merged/data/benchmark_metadata.csv"
  rank_files <- list.files("results/03_benchmark/merged/02_mean_rank/",
                            recursive = TRUE, full.names = T)
  gene_citations <- "resources/protein_citations.csv"
  overview_meta <- "results/manuscript_figures/figure_3/overview_kin.pdf"
  auroc_plot <- "results/manuscript_figures/figure_3/auroc_res.pdf"
  rank_plot <- "results/manuscript_figures/figure_3/mean_rank.pdf"
  kin_file <- "results/manuscript_figures/figure_3/kinase_GSknown.csv"
  heat_plot <- "results/manuscript_figures/figure_3/median_auroc.pdf"
  medRank_plot <- "results/manuscript_figures/figure_3/median_rank.pdf"
  cor_plot <- "results/manuscript_figures/figure_3/study_bias.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)
library(ComplexHeatmap)
library(circlize)

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
    add_column(net = str_split(file, "/")[[1]][5])
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

## AUROC Median heatmap
med_mat <- bench_df %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "score") %>%
  column_to_rownames("net")

col_fun = colorRamp2(c(min(med_mat, na.rm = T), max(med_mat, na.rm = T)), c("white", "deeppink4"))
column_ha = HeatmapAnnotation(mean_method = colMeans(med_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(med_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(heat_plot, width = 8, height = 4)
Heatmap(med_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", med_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

## Mean rank ---------------------------
priors <- map_chr(str_split(bench_files, "/"), 5) %>% unique()
rank_df <- map_dfr(rank_files, function(x) read_csv(x, col_types = cols()))
rank_df <- rank_df %>%
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

## scaled Median rank heatmap
medRank_mat <- rank_df %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(prior, method) %>%
  summarise(scaled_rank = mean(scaled_rank)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "scaled_rank") %>%
  column_to_rownames("prior")

col_fun = colorRamp2(c(min(medRank_mat), max(medRank_mat)), c("deepskyblue4", "white"))
column_ha = HeatmapAnnotation(mean_method = colMeans(medRank_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(medRank_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(medRank_plot, width = 8, height = 4)
Heatmap(medRank_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", medRank_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

## Mean rank per kinase ---------------------------
kin_gsknown <- rank_df %>%
  filter(prior == "phosphositeplus") %>%
  filter(!method == "number_of_targets") %>%
  group_by(targets) %>%
  summarise(mean_rank = mean(rank, na.rm = T)) %>%
  filter(!is.na(mean_rank)) %>%
  arrange(mean_rank) %>%
  mutate(mean_rank = round(mean_rank)) %>%
  dplyr::select(targets, mean_rank)

write_csv(kin_gsknown, kin_file)

## Compare to number of citations
overview_citations <- read_csv(gene_citations, col_types = cols())

kin_citations <- kin_gsknown %>%
  left_join(overview_citations %>% rename("targets" = Symbol), by = "targets")

cor_p <- ggplot(kin_citations, aes(x = mean_rank, y = log(citations))) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.6, fullrange = TRUE, color = "steelblue3") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.key.size = unit(0.4, 'cm'))  +
  geom_text(aes(label = targets), check_overlap = TRUE, size = 2.2, nudge_x = 0.1, nudge_y = 0.1)+
  ylab("log (# pubmed IDs)") +
  xlab("mean rank") +
  ggtitle(paste0("Pearson correlation: ", round(cor(kin_citations$mean_rank, log(kin_citations$citations)), digits = 2)))

cor.test(kin_citations$mean_rank, log(kin_citations$citations))

pdf(cor_plot, width = 4, height = 2.7)
cor_p
dev.off()


