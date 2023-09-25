#' Test overlap of priors

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$prior_files
  full_jaccard <- snakemake@output$full_jaccard
  kin_jaccard <- snakemake@output$kin_jaccard
}else{
  input_file <- list.files("results/hernandez/prior", full.names = T, pattern = ".tsv")
  full_jaccard <- "results/hernandez/overview_priors/overlap_priors.pdf"
  kin_jaccard <- "results/hernandez/overview_priors/overlap_priors_perKin.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ComplexHeatmap)

## Load data ---------------------------
priors <- map(input_file, read_tsv)
names(priors) <- str_remove(str_remove(input_file, "results/hernandez/prior/"), ".tsv")

priors <- priors[!str_detect(names(priors), "_")]
priors <- priors[!names(priors) == "jhonson"]

## Extract links ---------------------------
links_priors <- map(priors, function(prior){
  links <- prior %>%
    mutate(target = str_remove(target, "\\|auto")) %>%
    filter(str_detect(target, "\\|")) %>%
    mutate(link = paste(source, target, sep = ";"))
  kin_msk <- links %>%
    group_by(source) %>%
    summarise(n = n()) %>%
    filter(n >= 5) %>%
    pull(source)
  links %>%
    filter(source %in% kin_msk) %>%
    pull(link)
})

dist <- unlist(lapply(combn(links_priors, 2, simplify = FALSE), function(x) {
  length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))

jacc_links <- cbind(t(combn(names(links_priors),2)), dist) %>%
  as.data.frame() %>%
  mutate(dist = as.numeric(dist))

jacc_links <- rbind(jacc_links,
                    jacc_links %>%
                     mutate(V1_tmp = V2) %>%
                     mutate(V2 = V1) %>%
                     mutate(V1 = V1_tmp) %>%
                     select(-V1_tmp),
                   data.frame(V1 = unique(c(jacc_links$V1, jacc_links$V2)),
                              V2 = unique(c(jacc_links$V1, jacc_links$V2)),
                              dist = 1))
jacc_links_wide <- jacc_links %>%
  pivot_wider(names_from = V2, values_from = dist) %>%
  column_to_rownames("V1")

jacc_links_wide <- jacc_links_wide[colnames(jacc_links_wide), colnames(jacc_links_wide)]

heatmap_fullJacc <- Heatmap(jacc_links_wide)


## Jacc per kinase ---------------------------
kinases <- map(priors, function(i) i %>% pull(source)) %>%
  unlist() %>%
  unique()

jacc_per_kin <- map_dfr(kinases, function(kin_i){
  if((match(kin_i, kinases)/50)%%1==0){
    print(paste0(match(kin_i, kinases), "/", length(kinases)))
  }
  links_kin <- map(priors, function(prior){
    links <- prior %>%
      mutate(target = str_remove(target, "\\|auto")) %>%
      filter(str_detect(target, "\\|")) %>%
      filter(source == kin_i) %>%
      mutate(link = paste(source, target, sep = ";"))

    if (nrow(links) >= 5){
      links %>%
        pull(link)
    } else {
      c(NA)
    }
  })

  dist <- unlist(lapply(combn(links_kin, 2, simplify = FALSE), function(x) {
    if (!any(is.na(x[[1]])) & !any(is.na(x[[2]]))){
      length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]]))
      } else {
      NA
        }
  }))
  cbind(t(combn(names(links_kin),2)), dist) %>%
    as.data.frame() %>%
    mutate(dist = as.numeric(dist)) %>%
    add_column(kinase = kin_i)
})

mean_jacc <- jacc_per_kin %>%
  group_by(V1, V2) %>%
  summarise(mean_jaccard = mean(dist, na.rm = T)) %>%
  ungroup()

mean_jacc <- rbind(mean_jacc,
      mean_jacc %>%
        mutate(V1_tmp = V2) %>%
        mutate(V2 = V1) %>%
        mutate(V1 = V1_tmp) %>%
        select(-V1_tmp),
      data.frame(V1 = unique(c(mean_jacc$V1, mean_jacc$V2)),
                 V2 = unique(c(mean_jacc$V1, mean_jacc$V2)),
                 mean_jaccard = 1))

mean_jac_wide <- mean_jacc %>%
  pivot_wider(names_from = V2, values_from = mean_jaccard) %>%
  column_to_rownames("V1")

mean_jac_wide <- mean_jac_wide[colnames(mean_jac_wide), colnames(mean_jac_wide)]

heatmap_kinJacc <- Heatmap(mean_jac_wide)

## Save plots ---------------------------
pdf(full_jaccard, height = 5, width = 5)
heatmap_fullJacc
dev.off()

pdf(kin_jaccard, height = 5, width = 5)
heatmap_kinJacc
dev.off()
