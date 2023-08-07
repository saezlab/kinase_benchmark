#' Get mean/median rank for each method

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  output_file <- snakemake@output$output
  meta_out <- snakemake@output$meta_out
}else{
  input_files <- list.files("results/decryptm/activity_scores", full.names = T)
  meta_file <- "results/decryptm/processed_data/meta_data.csv"
  output_file <- "results/decryptm/benchmark/mean_rank.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

## Load  meta ---------------------------
obs <- read_csv(meta_file) %>%
  dplyr::select(drug, dose, dataset, cells, sample, Target)

## Get gene names of targets in meta ---------------------------
all_targets <- str_split(obs$Target, ";") %>% unlist()
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res <- getBM(attributes = c('uniprot_gn_id',
                            'external_gene_name'),
             values = all_targets,
             mart = mart)

obs_targets <- map_dfr(1:nrow(obs), function(i){
  uniprotIDs <- str_split(obs[i,]$Target, ";") %>% unlist()

  target_gene <- res %>%
    filter(uniprot_gn_id %in% uniprotIDs) %>%
    pull(external_gene_name)

  obs[i,] %>% add_column(perturb = paste(target_gene, collapse = ";"))
})

# filter out experiments with unknown target (e.g. several members)
target_df <- obs_targets %>%
  dplyr::filter(!perturb == "") %>%
  add_column(sign = 1)

## Load  activity scores ---------------------------
ranks <- map_dfr(input_files, function(input_file){
  net <- str_remove(str_split(str_split(input_file, "/")[[1]][4], "-")[[1]][2], ".rds")
  act_scores <- readRDS(input_file)
    map_dfr(names(act_scores), function(meth){
      method_act <- act_scores[[meth]]

      method_act_long <- method_act %>%
        rownames_to_column("kinase") %>%
        pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
        filter(!is.na(score)) %>%
        arrange(score)

      map_dfr(unique(method_act_long$sample), function(exp){
        act_df <- method_act_long %>%
          dplyr::filter(sample == exp)

        if (exp %in% target_df$sample){
          targets <- target_df %>%
            dplyr::filter(sample == exp) %>%
            pull(perturb) %>%
            str_split(";") %>%
            unlist()


          rank <- map_dbl(targets, function(target){
            position <- which(act_df$kinase %in% target)
            if (length(position) == 0){
              position <- NA
            }
            position
          })

          data.frame(sample = exp,
                     method = meth,
                     prior = net,
                     targets = targets,
                     rank = rank,
                     kinases = nrow(act_df)) %>%
            mutate(scaled_rank = rank/kinases)
        }

      })

    })
})

mean_rank_df <- ranks %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(method, prior) %>%
  summarise(mean_rank = round(mean(rank), digits = 1), mean_scaled_rank = round(mean(scaled_rank), digits = 2)) %>%
  arrange(mean_scaled_rank, mean_rank) %>%
  filter(!is.na(prior))

write_csv(mean_rank_df, output_file)

