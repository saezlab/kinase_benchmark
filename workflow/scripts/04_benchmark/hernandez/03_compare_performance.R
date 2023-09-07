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
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = order_m)


auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=method)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(hjust=0.5)) +
  ylab("AUROC") +
  xlab("")

auprc_p <- bench_df %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=method)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(hjust=0.5)) +
  ylab("AUPRC") +
  xlab("")

pdf(auroc_plot, width = 10, height = 6)
auroc_p
dev.off()

pdf(auprc_plot, width = 10, height = 6)
auprc_p
dev.off()
