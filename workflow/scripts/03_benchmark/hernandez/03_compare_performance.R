#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  auroc_plot <- snakemake@output$auroc
  auprc_plot <- snakemake@output$auprc
}else{
  bench_files <- list.files("results/hernandez/benchmark_res",
                           pattern = "bench", recursive = TRUE, full.names = T)
  auroc_plot <- "results/hernandez/benchmark_res/plots/AUROC.pdf"
  auprc_plot <- "results/hernandez/benchmark_res/plots/AUPRC.pdf"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
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


auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=1, lwd=0.6) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.5)) +
  ylab("AUROC") +
  xlab("")

auprc_p <- bench_df %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=3, lwd=1) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 0.5)) +
  ylab("AUPRC") +
  xlab("")

pdf(auroc_plot, width = 12, height = 5)
auroc_p
dev.off()

pdf(auprc_plot, width = 12, height = 5)
auprc_p
dev.off()
