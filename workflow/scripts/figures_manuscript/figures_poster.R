#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  auroc_plot <- snakemake@output$auroc
  auprc_plot <- snakemake@output$auprc
}else{
  GPSppsp_files <- list.files("results/hernandez/benchmark_res/GPSppsp/scaled",
                            pattern = "bench", recursive = TRUE, full.names = T)
  all_files <- list.files("results/hernandez/benchmark_res",
                              pattern = "bench", recursive = TRUE, full.names = T)
  auroc_plot <- "results/manuscript_figures/poster/AUROC_GPSppsp.png"
  auroc_plot_so <- "results/manuscript_figures/poster/AUROC_GPSppsp_same_order.png"
  auroc_plot_combined <- "results/manuscript_figures/poster/AUROC_predTargets.png"
  auroc_plot_combined_so <- "results/manuscript_figures/poster/AUROC_predTargets_same_order.png"
  GPSnetworkin_plot_combined <- "results/manuscript_figures/poster/AUROC_GPSnetworkin.png"
  GPSnetworkin_plot_combined_so <- "results/manuscript_figures/poster/AUROC_GPSnetworkin_same_order.png"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
bench_list <- map(GPSppsp_files, function(file){
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


auroc_gpsppsp <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill = method)) +
  stat_boxplot(geom ='errorbar', width=0.4, lwd = 0.5, linetype = "longdash")+
  geom_boxplot(outlier.size=3, lwd=0.7, outlier.shape = 1, outlier.fill = "white") +
  scale_fill_manual(values=c( "blue", "yellow", "red", "orange",
                              "purple", "green", "lightblue", "pink",
                              "blue", "yellow", "red", "orange",
                              "purple", "green", "lightblue", "pink",
                              "blue", "yellow", "red")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("AUROC") +
  xlab("") +
  ylim(0.5, 0.9)

png(auroc_plot, width = 1000, height = 300)
auroc_gpsppsp
dev.off()
pdf(str_replace(auroc_plot, ".png", ".pdf"), width = 10, height = 4)
auroc_gpsppsp
dev.off()

bench_df$method <- factor(bench_df$method, levels = c("RoKAI_z", "fgsea", "mean", "wsum",
                                                      "ptmsea", "norm_fgsea", "viper", "norm_wmean",
                                                      "KSEA_z", "ulm", "Wilcox", "PC1",
                                                      "rokai_lm", "UQ", "KS", "median",
                                                      "mlm", "KARP", "number_of_targets"))

auroc_gpsppsp <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill = method)) +
  stat_boxplot(geom ='errorbar', width=0.4, lwd = 0.7, linetype = "longdash")+
  geom_boxplot(outlier.size=3, lwd=0.5, outlier.shape = 1, outlier.fill = "white") +
  scale_fill_manual(values=c( "blue", "yellow", "red", "orange",
                              "purple", "green", "lightblue", "pink",
                              "blue", "yellow", "red", "orange",
                              "purple", "green", "lightblue", "pink",
                              "blue", "yellow", "red")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("AUROC") +
  xlab("") +
  ylim(0.5, 0.9)

png(auroc_plot_so, width = 1000, height = 370)
auroc_gpsppsp
dev.off()

pdf(str_replace(auroc_plot_so, ".png", ".pdf"), width = 10, height = 4)
auroc_gpsppsp
dev.off()


## Known + predicted targets
files <- all_files[str_detect(all_files, "GPSppsp")]
files <- files[str_detect(files, "scaled_subset")]

bench_list <- map(files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][4])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)%>%
  filter(method %in% c("RoKAI_z", "mean", "wsum", "ptmsea", "fgsea", "viper", "KSEA_z", "number_of_targets"))


bench_df$net <- factor(bench_df$net, levels = c("GPSppsp", "GPSppsp_networkin", "GPSppsp_iKiPdb"))
bench_df$method[bench_df$method == "number_of_targets"] <- "control"
bench_df$method[bench_df$method == "KSEA_z"] <- "KSEA"

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = order_method)


p1 <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill = net)) +
  geom_boxplot(outlier.size=3, lwd=0.7, outlier.shape = 1, outlier.fill = "white") +
  scale_fill_manual(values=c( "#F8766D", "#00BA38", "#619CFF")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("AUROC") +
  xlab("") +
  ylim(0.5, 0.9)

bench_df$method <- factor(bench_df$method, levels = c("RoKAI_z", "mean", "wsum", "ptmsea", "fgsea", "KSEA", "viper", "control"))


p2 <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill = net)) +
  geom_boxplot(outlier.size=3, lwd=0.7, outlier.shape = 1, outlier.fill = "white") +
  scale_fill_manual(values=c( "#F8766D", "#00BA38", "#619CFF")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("AUROC") +
  xlab("") +
  ylim(0.5, 0.9)

png(auroc_plot_combined, width = 600, height = 300)
p1
dev.off()

pdf(str_replace(auroc_plot_combined, ".png", ".pdf"), width = 6, height = 3)
p1
dev.off()

png(auroc_plot_combined_so, width = 600, height = 300)
p2
dev.off()

pdf(str_replace(auroc_plot_combined_so, ".png", ".pdf"), width = 6, height = 3)
p2
dev.off()


## GPSppsp_Network
files <- all_files[str_detect(all_files, "GPSppsp_networkin")]
files <- files[str_detect(files, "scaled")]

bench_list <- map(files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][4])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)%>%
  filter(method %in% c("RoKAI_z", "mean", "wsum", "ptmsea", "fgsea", "viper", "KSEA_z", "number_of_targets"))


bench_df$method[bench_df$method == "number_of_targets"] <- "control"
bench_df$method[bench_df$method == "KSEA_z"] <- "KSEA"

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = order_method)


p1 <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill = method)) +
  geom_boxplot(outlier.size=3, lwd=0.7, outlier.shape = 1, outlier.fill = "white") +
  scale_fill_manual(values=c( "blue", "yellow", "red", "orange",
                              "purple", "green", "lightblue", "pink")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("AUROC") +
  xlab("") +
  ylim(0.5, 0.9)

bench_df$method <- factor(bench_df$method, levels = c("RoKAI_z", "mean", "wsum", "ptmsea", "fgsea", "KSEA", "viper", "control"))

p2 <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill = method)) +
  geom_boxplot(outlier.size=3, lwd=0.7, outlier.shape = 1, outlier.fill = "white") +
  scale_fill_manual(values=c( "blue", "yellow", "red", "orange",
                              "purple", "green", "lightblue", "pink")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("AUROC") +
  xlab("") +
  ylim(0.5, 0.9)


png(GPSnetworkin_plot_combined, width = 600, height = 300)
p1
dev.off()

pdf(str_replace(GPSnetworkin_plot_combined, ".png", ".pdf"), width = 6, height = 3)
p1
dev.off()

png(GPSnetworkin_plot_combined_so, width = 600, height = 300)
p2
dev.off()

pdf(str_replace(GPSnetworkin_plot_combined_so, ".png", ".pdf"), width = 6, height = 3)
p2
dev.off()
