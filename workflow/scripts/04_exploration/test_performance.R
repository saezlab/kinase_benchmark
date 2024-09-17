#' Get mean/median rank for each method

## Snakemake ---------------------------
if(exists("snakemake")){
  input_files <- snakemake@input$scores
  out_plot <- snakemake@output$out
}else{
  input_files <- "results/03_benchmark/hernandez/02_mean_rank/GPS/mean-GPS.csv"
  out_plot <- "results/03_benchmark/hernandez/02_mean_rank/plots/mean-GPS_targets.pdf"
}

## Libraries ---------------------------
library(tidyverse)

## Load scores and meta ---------------------------
rank_df <- read_csv(input_files, col_types = cols())

net <- str_split(str_remove(str_split(input_files, "/")[[1]][6], ".csv"), "-")[[1]][2]
meth <- str_split(str_remove(str_split(input_files, "/")[[1]][6], ".csv"), "-")[[1]][1]

## Get rank ---------------------------
rank_df <- rank_df %>%
  group_by(targets, target_size, class) %>%
  summarise(mean_scaled_rank = mean(scaled_rank, na.rm = TRUE),
            median_scaled_rank = median(scaled_rank, na.rm = TRUE))

p1 <- ggplot(rank_df, aes(x = target_size, y = mean_scaled_rank)) +
  geom_boxplot() +
  labs(title = paste0(net, "-", meth, " target size"),
       x = "Target size",
       y = "Scaled rank") +
  theme_minimal()

p2 <- ggplot(rank_df, aes(x = class, y = mean_scaled_rank)) +
  geom_boxplot() +
  labs(title = paste0(net, "-", meth, " kinase class"),
       x = "Target size",
       y = "Scaled rank") +
  theme_minimal()

table(rank_df$target_size)
table(rank_df$class)
## Get rank ---------------------------
pdf(out_plot, width = 5, height = 5)
p1
p2
dev.off()
