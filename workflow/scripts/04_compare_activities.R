if(exists("snakemake")){
  act_files <- snakemake@input$act_files
  plot_dest <- snakemake@output$pdf
}else{
  act_files <- list.files("results/activity_scores", pattern = "rds", recursive = T, full.names = T)
  plot_dest <- "results/comparison/plots/Pearson_heatmap.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ComplexHeatmap)

## Load and merge scores ---------------------------
act_list <- map(act_files, readRDS)
names(act_list) <- str_remove(map_chr(str_split(act_files, "/"), 4), ".rds")
act_df <- map_dfr(names(act_list), function(act_i){
  act <- act_list[[act_i]]
  act_df <- map_dfr(names(act), function(method){
    act_method <- act[[method]]
    act_method %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      add_column(method = method)
  })
  act_df %>%
    add_column(prior = str_split(act_i, "_")[[1]][2]) %>%
    add_column(cancer = str_split(act_i, "_")[[1]][1])
})

act_df <- act_df %>%
  group_by(cancer)

act_cancer <- act_df %>%
  group_split()

names(act_cancer) <- group_keys(act_df)$cancer

correlate_scores <- function(df, long = T){
  df_wide <- df %>%
    mutate(id_method = paste(method, prior, sep = ":")) %>%
    mutate(id_sample = paste(kinase, sample, sep = ":")) %>%
    dplyr::select(id_method, id_sample, score) %>%
    pivot_wider(names_from = id_method, values_from = score) %>%
    column_to_rownames("id_sample")
  df_wide[is.infinite(as.matrix(df_wide))] <- NA

  cor_df <- cor(df_wide, method = "pearson", use = "pairwise.complete.obs")

  if (long){
    df_long <- cor_df %>%
      as.data.frame() %>%
      rownames_to_column("id") %>%
      pivot_longer(!id, names_to = "id2", values_to = "Pearson") %>%
      mutate(cancer = unique(df$cancer))

    return(df_long)
  } else {
    return(cor_df)
  }
}

cor_scores <- map_dfr(act_cancer, correlate_scores)

meanPearson <- cor_scores %>%
  group_by(id, id2) %>%
  summarise(mean_Pearson = mean(Pearson)) %>%
  pivot_wider(names_from = "id2", values_from = "mean_Pearson") %>%
  column_to_rownames("id")

column_ha <- HeatmapAnnotation(method = map_chr(str_split(colnames(meanPearson), ":"), 1),
                              prior = map_chr(str_split(colnames(meanPearson), ":"), 2))

# note how we set the width of this empty annotation
cor_plot <- Heatmap(meanPearson,
        name = "Pearson",
        show_row_names = F,
        show_column_names = F,
        top_annotation = column_ha)

pdf(file=plot_dest, height = 10, width = 10)
draw(cor_plot)
dev.off()
