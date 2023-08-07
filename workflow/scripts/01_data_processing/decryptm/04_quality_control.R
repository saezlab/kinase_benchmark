#'

## Snakemake ---------------------------
if(exists("snakemake")){
  EC50_file <- snakemake@input$EC50_file
  output_plots <- snakemake@output$output_plots
}else{
  EC50_file <- "results/decryptm/processed_data/R2_pEC50.csv"
  output_plots <- "results/decryptm/processed_data/figures/QC_R2_pEC50.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(pheatmap)
library(ggVennDiagram)
library(ggfortify)
library(ggplot2)
library(grid)

## Load data ---------------------------
EC50 <- read_csv(EC50_file) %>%
  column_to_rownames("pps_id")

meta <- read_csv("results/decryptm/processed_data/meta_data.csv") %>%
  rowwise %>%
  mutate(tmp_id = paste(c(drug, dose, dose_scale, time, cells), collapse = ":"))

## QC ---------------------------
## number of phosphorylation sites with assigned EC50 values
n_ec50 <- colSums(!is.na(EC50)) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  arrange(desc(.))
n_ec50$sample <- factor(n_ec50$sample, levels = unique(n_ec50$sample))

p1 <- ggplot(n_ec50, aes(x = sample, y = .)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of phosphorylation sites\nwith assigned EC50 values")

## EC50 value distribution per sample
ed50_long <- EC50 %>%
  rownames_to_column("pps") %>%
  pivot_longer(!pps, names_to = "sample", values_to = "EC50")
ed50_long$sample <- factor(ed50_long$sample, levels = unique(n_ec50$sample))
p2 <- ggplot(ed50_long, aes(x = sample, y = abs(EC50))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("-log(EC50)")

p2.2 <- ggplot(ed50_long, aes(x = sample, y = EC50)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("-log(EC50) * sign(effect size)")

# Phosphosite NA distribution
na_dist <- is.na(EC50) %>%
  rowSums() %>%
  as.data.frame() %>%
  add_column(n = ncol(EC50)) %>%
  mutate(perc = 1 - ./n)

p2.3 <- ggplot(na_dist, aes(x=perc)) +
  geom_histogram(bins = 50) +
  xlab("PPS measured across sampels [%]")

# correlation matrix
cor_m <- cor(EC50, use = "pairwise.complete.obs")

p3 <- pheatmap(cor_m, show_rownames = F)

# Venn Diagrams for replicates
replicates <- meta$tmp_id[duplicated(meta$tmp_id)] %>% unique()

p4_list <- map(replicates, function(rep_i){
  sample_rep <- meta %>% filter(tmp_id == rep_i) %>% pull(sample)
  EC50_rep <- EC50[,sample_rep]

  list_pps <- map(colnames(EC50_rep), function(col_i){
    rownames(EC50_rep)[!is.na(EC50_rep[,col_i])]
  })
  names(list_pps) <- colnames(EC50_rep)

  ggVennDiagram(list_pps)
})

# PCA
EC50_filtered <- EC50[rowSums(is.na(EC50)) < 0.6*ncol(EC50),]
EC50_filtered[is.na(EC50_filtered)] <- 0


pca_res <- prcomp(t(EC50_filtered), scale. = TRUE)

p5 <- ggplot2::autoplot(pca_res, data = cbind(t(EC50_filtered), data.frame(sample = colnames(EC50))), colour = "sample")

legend <- cowplot::get_legend(p5)

grid.newpage()
grid.draw(legend)

# histograms per sample
p6_list <- map(colnames(EC50), function(col_i){
  # Some index conversion as otherwise there is a problem with one column name
  # Error in parse(text = paste_line(x)) : <text>:1:11: unexpected input 1: ddPTM_PC-9_
  idx <- which(colnames(EC50) == (col_i))
  tmp <- EC50
  colnames(tmp)[idx] <- "EC50_values"

  #Plot histogram
  ggplot(tmp, aes(x = EC50_values)) +
    geom_histogram(bins = 50) +
    xlab("-log(EC50)") +
    ggtitle(col_i)
})


## Save plots ---------------------------
pdf(file = output_plots, height = 15, width = 15)
p1
p2
p2.2
p2.3
grid.newpage()
p3
p5 + theme(legend.position = "none")
grid.newpage()
grid.draw(legend)
p4_list
p6_list
dev.off()

