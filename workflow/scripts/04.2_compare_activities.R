if(exists("snakemake")){
  act_files <- snakemake@input$act_files
  jaccard_i <-  snakemake@params$jaccard_i
  plot_spearman <- snakemake@output$plotSpearman
  plot_pearson <-  snakemake@output$plotPearson
  plot_jaccard_up <-  snakemake@output$plotJaccardUp
  plot_jaccard_down <-  snakemake@output$plotJaccardDown
  dist_csv <- snakemake@output$dist_csv
}else{
  act_files <- list.files("results/final_scores", pattern = "rds", recursive = T, full.names = T)
  plot_spearman <- "results/comparison/plots/spearman_heatmap.pdf"
  plot_pearson <- "results/comparison/plots/pearson_heatmap.pdf"
  plot_jaccard_up <- "results/comparison/plots/jaccard_up_heatmap.pdf"
  plot_jaccard_down <- "results/comparison/plots/jaccard_down_heatmap.pdf"
  jaccard_i <- 10
  dist_csv <- "results/comparison/mean_distance.csv"
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
    act_method <- act_method %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      add_column(method = method)

    if(method == "mlm"){
      act_method <- map_dfr(1:nrow(act_method), function(row_i){
        df <- act_method[row_i,]
        if (!is.null(unlist(df$score))){
          df %>%
            mutate(score = unlist(df$score)[1])
        } else if (is.null(unlist(df$score))){
          df %>%
            mutate(score = NA)
        }
      })
    }
    act_method
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

## Score comparison ---------------------------
correlate_scores <- function(df, method = "spearman", long = T){
  df_wide <- df %>%
    mutate(id_method = paste(method, prior, sep = ":")) %>%
    mutate(id_sample = paste(kinase, sample, sep = ":")) %>%
    dplyr::select(id_method, id_sample, score) %>%
    pivot_wider(names_from = id_method, values_from = score) %>%
    column_to_rownames("id_sample")
  df_wide[is.infinite(as.matrix(df_wide))] <- NA

  cor_df <- cor(df_wide, method = method, use = "complete.obs")

  if (long){
    df_long <- cor_df %>%
      as.data.frame() %>%
      rownames_to_column("id") %>%
      pivot_longer(!id, names_to = "id2", values_to = "statistic") %>%
      mutate(cancer = unique(df$cancer))

    return(df_long)
  } else {
    return(cor_df)
  }
}

jaccard_index <- function(df_list, scores = "up", n = 10){
  top <- map_dfr(df_list, function(df){
    df$score[is.infinite(df$score)] <- NA
    df <- df %>%
      na.omit() %>%
      arrange(desc(score))

    df_sample_list <- df %>%
      group_by(sample) %>%
      group_split()

    jaccard_sample <- map_dfr(df_sample_list, function(df_sample){
      if (scores == "up") {
        kin_df <- df_sample %>%
          mutate(id = paste(prior, method, sep = ":")) %>%
          group_by(id) %>%
          arrange(desc(score)) %>%
          slice(1:n)
        } else if (scores == "down"){
          kin_df <- df_sample %>%
            mutate(id = paste(prior, method, sep = ":")) %>%
            group_by(id) %>%
            arrange(score) %>%
            slice(1:n)
        }

      kin_list <- kin_df %>%
        group_split()
      names(kin_list) <- group_keys(kin_df)$id

      kin <- map(kin_list, function(df){
        df %>% pull(kinase)
        })

      dist <- unlist(lapply(combn(kin, 2, simplify = FALSE), function(x) {
        length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))
      jaccard_df <- cbind(data.frame(t(combn(names(kin), 2))), dist)

      # duplicate to have all values
      full_jaccard <- rbind(jaccard_df,
                            data.frame(X1 = jaccard_df$X2,
                                       X2 = jaccard_df$X1,
                                       dist = jaccard_df$dist),
                            data.frame(X1 = unique(jaccard_df$X1),
                                       X2 = unique(jaccard_df$X1),
                                       dist = 1))

      full_jaccard %>%
        mutate(sample = unique(kin_df$sample)) %>%
        mutate(cancer = unique(kin_df$cancer)) %>%
        rename("jaccard" = dist) %>%
        arrange(desc(jaccard))
      })

    jaccard_sample
    })

  top_summarised <- top %>%
    group_by(X1, X2) %>%
    summarise(jaccard_mean = mean(jaccard), .groups = "keep") %>%
    arrange(desc(jaccard_mean))

  return(top_summarised)
}

## Spearman correlation ---------
cor_scores_spearman <- map_dfr(act_cancer, correlate_scores)

meanSpearman <- cor_scores_spearman %>%
  group_by(id, id2) %>%
  summarise(mean_Spearman = mean(statistic)) %>%
  pivot_wider(names_from = "id2", values_from = "mean_Spearman") %>%
  column_to_rownames("id")

column_ha <- HeatmapAnnotation(method = map_chr(str_split(colnames(meanSpearman), ":"), 1),
                              prior = map_chr(str_split(colnames(meanSpearman), ":"), 2))

cor_plot_spearman <- Heatmap(as.matrix(meanSpearman),
        name = "Spearman",
        show_row_names = F,
        show_column_names = F,
        top_annotation = column_ha)

pdf(file=plot_spearman, height = 10, width = 10)
draw(cor_plot_spearman)
dev.off()

## Pearson correlation --------
cor_scores_pearson <- map_dfr(act_cancer, function(df){correlate_scores(df = df, method = "pearson")})

meanPearson <- cor_scores_pearson %>%
  group_by(id, id2) %>%
  summarise(mean_Pearson = mean(statistic)) %>%
  pivot_wider(names_from = "id2", values_from = "mean_Pearson") %>%
  column_to_rownames("id")

column_ha <- HeatmapAnnotation(method = map_chr(str_split(colnames(meanPearson), ":"), 1),
                               prior = map_chr(str_split(colnames(meanPearson), ":"), 2))

cor_plot_pearson <- Heatmap(as.matrix(meanPearson),
                    name = "Pearson",
                    show_row_names = F,
                    show_column_names = F,
                    top_annotation = column_ha)

pdf(file=plot_pearson, height = 10, width = 10)
draw(cor_plot_pearson)
dev.off()

## Jaccard ---------------------------
jaccard_i <- as.numeric(jaccard_i)
jaccard_up <- jaccard_index(act_cancer, n = jaccard_i)
jaccard_down <- jaccard_index(act_cancer, scores = "down", n = jaccard_i)

jaccard_up <- jaccard_up %>%
  pivot_wider(names_from = "X2", values_from = "jaccard_mean") %>%
  column_to_rownames("X1")

column_ha <- HeatmapAnnotation(method = map_chr(str_split(colnames(jaccard_up), ":"), 2),
                               prior = map_chr(str_split(colnames(jaccard_up), ":"), 1))

jaccard_plot_up <- Heatmap(as.matrix(jaccard_up),
                            name = "Jaccard",
                            show_row_names = F,
                            show_column_names = F,
                            top_annotation = column_ha)

pdf(file=plot_jaccard_up, height = 10, width = 10)
draw(jaccard_plot_up)
dev.off()

jaccard_down <- jaccard_down %>%
  pivot_wider(names_from = "X2", values_from = "jaccard_mean") %>%
  column_to_rownames("X1")

column_ha <- HeatmapAnnotation(method = map_chr(str_split(colnames(jaccard_down), ":"), 1),
                               prior = map_chr(str_split(colnames(jaccard_down), ":"), 2))

jaccard_plot_down <- Heatmap(as.matrix(jaccard_down),
                           name = "Jaccard",
                           show_row_names = F,
                           show_column_names = F,
                           top_annotation = column_ha)

pdf(file=plot_jaccard_down, height = 10, width = 10)
draw(jaccard_plot_down)
dev.off()

## Calculate effect of method and prior knowledge on result ---------------------------
distSpearman <- dist(meanSpearman) %>%
  as.matrix()
distSpearman[upper.tri(meanSpearman)] <- NA
distSpearman <- distSpearman %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "id2", values_to = "dist") %>%
  dplyr::filter(!id == id2) %>%
  dplyr::filter(!is.na(dist))

method_PK <- unique(unlist(str_split(distSpearman$id, ":")))

effect_method_PKN <- map_dfr(method_PK, function(x){
  x_raw <- x
  if (x == "INKA"){
    x <- "INKA:" #to avoid confusion with INKA_kinase_centric
  }
  if (x == "wmean"){
    x <- "\\bwmean" #to avoid confusion with norm_wmean
  }
  distance_within <- distSpearman %>%
    dplyr::filter(str_detect(id, x) & str_detect(id2, x)) %>%
    pull(dist) %>%
    mean()

  distance_withothers <- distSpearman %>%
    dplyr::filter(str_detect(id, x) | str_detect(id2, x)) %>%
    dplyr::filter(!(str_detect(id, x) & str_detect(id2, x))) %>%
    pull(dist) %>%
    mean()

  data.frame(method_PKN = x_raw,
             dist_within = distance_within,
             dist_to_others = distance_withothers,
             dist_ratio = distance_within/distance_withothers)
})

effect_method_PKN <- effect_method_PKN %>%
  arrange(dist_ratio)

write_csv(effect_method_PKN, file = dist_csv)
