#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  output_file <- snakemake@output$output
  meta_out <- snakemake@output$meta_out
}else{
  input_file <- "results/decryptm/activity_scores/R2_pEC50-GPS.rds"
  meta_file <- "results/decryptm/processed_data/meta_data.csv"
  meta_out <- "results/decryptm/benchmark_scores/obs_INKA-R2_pEC50-GPS.csv"
  output_file <- "results/decryptm/benchmark_scores/INKA-R2_pEC50-GPS.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

method_tmp <- str_remove(str_split(output_file, "/")[[1]][4], ".csv")
input <- paste0("-",str_remove(str_split(input_file, "/")[[1]][4], ".rds"))
method <- gsub(input, "", method_tmp)
print(method)

## Load scores and meta ---------------------------
act_scores <- readRDS(input_file)

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

## Get scores for each method ---------------------------
# change direction to perturbation
df <- t(act_scores[[method]]) * -1
df <- df %>%
  as.data.frame() %>%
  rownames_to_column("experiment")

df <- df[df$experiment %in% target_df$sample,]

# filter out experiments where no activity was infurred for perturbed kinases
df_filtered <- map_dfr(1:nrow(df), function(i){
  tmp <- df[i,]
  targets <- target_df %>%
    filter(sample %in% tmp$experiment) %>%
    pull(perturb) %>%
    str_split(";") %>%
    unlist

  sum_act <- sum(tmp[,colnames(tmp) %in% targets])

  if(!is.na(sum_act) & !sum_act == 0){
    df[i,]
  } else {
    df[i,] %>% mutate(experiment = "remove")
  }
})

df_filtered <- df_filtered %>%
  filter(!experiment == "remove")

write_csv(df_filtered, output_file)

# filter out meta to fit to experiments in matrix
target_df <- target_df %>%
  filter(sample %in% df_filtered$experiment)

write_csv(target_df, meta_out)

