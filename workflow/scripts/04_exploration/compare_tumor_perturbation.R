#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  tumor_files <- snakemake@input$tumor
  performance_plot <- snakemake@output$plot
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res_subset/subset",
                            pattern = "bench", recursive = TRUE, full.names = T)
  tumor_files <- "data/misc/known_roc_filt.Rds"                  
  performance_plot <- "results/04_exploration/all/performance_subset.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## Load AUROC ---------------------------
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
df_perturb <- bind_rows(bench_list) %>%
  dplyr::filter(net == "phosphositeplus") %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(method, score) %>%
  add_column(benchmark = "perturbation-based")

known_roc <- readRDS(tumor_files)
roc <- lapply(known_roc$ROC_results, "[[", "sample_AUROCs")

df_tumor <- do.call(rbind, lapply(names(roc), function(method) {
  data.frame(method = method, score = roc[[method]], benchmark = "tumor-based")
}))

bench_df <- rbind(df_perturb, df_tumor)

mean_auroc <- bench_df %>%
  group_by(method, benchmark) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc))  # Sorting methods by mean AUROC

# Update df to ensure methods are ordered based on mean AUROC
bench_df <- bench_df %>%
  mutate(method = factor(method, levels = mean_auroc$method))

# Create the boxplot with ggplot2
lines <- mean_auroc$mean_auroc
names(lines) <- mean_auroc$method
auroc_p <- ggplot(bench_df, aes(x = method, y = score, fill = benchmark))  +
  geom_boxplot(outlier.size = 1, width = 1) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +  # Muted scientific color palette
  theme_minimal() +
  facet_grid(~method, scale='free_x') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.text.x = element_blank()) + 
    geom_hline(aes(yintercept = lines[method]), color = "#9E2710")


# Test scatter plot
mean_auroc_df <- bench_df %>%
  group_by(method, benchmark) %>%
  summarize(mean_auroc = mean(score)) %>%
  pivot_wider(names_from = "benchmark", values_from = "mean_auroc")

scatter_p <- ggplot(mean_auroc_df, aes(x=`perturbation-based`, y = `tumor-based`, color = method)) +
  geom_point() +
  theme_minimal()
    
pdf(performance_plot, width = 8, height = 4)
auroc_p 
scatter_p
dev.off()
