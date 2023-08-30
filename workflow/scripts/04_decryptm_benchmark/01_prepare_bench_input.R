#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  meta_file <- snakemake@input$meta
  processed_meta <- snakemake@input$kinome
  perturb <- snakemake@params$perturb
  output_file <- snakemake@output$output
  meta_out <- snakemake@output$meta_out
}else{
  input_file <- "results/decryptm/activity_scores/R2_pEC50-GPS.rds"
  meta_file <- "results/decryptm/processed_data/meta_data.csv"
  processed_meta <- "data/decryptm/decryptm_processed_targets.csv"
  perturb <- "kinomebeads"
  meta_out <- "results/decryptm/benchmark_scores/obs_INKA-R2_pEC50-GPS.csv"
  output_file <- "results/decryptm/benchmark_scores/INKA-R2_pEC50-GPS.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

method_tmp <- str_remove(str_split(output_file, "/")[[1]][4], ".csv")
input <- paste0("-",str_remove(str_split(input_file, "/")[[1]][4], ".rds"))
method <- gsub(input, "", method_tmp)

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

  obs[i,] %>% add_column(perturb = paste(unique(target_gene), collapse = ";"))
})

processed_targets <- read.csv(processed_meta)
targets <- processed_targets %>%
  mutate(drug = recode(drug,
                       "MK-2206" = "MK2206",
                       "Mirdametinib" = "PD325901",
                       "AZD-8055" = "AZD8055",
                       "Nitedanib" = "Nintedanib",
                       "Tideglusib" = "Tideglusib",
                       "Staurosporine" = "Staursporin",
                       "AZD-4547" = "AZD4547")) %>%
  group_by(drug) %>%
  mutate(targets = paste0(unique(gene_name), collapse = ";")) %>%
  dplyr::select(drug, targets) %>%
  distinct()

obs_targets <- left_join(obs_targets, targets, by = "drug")

# Combine targets for drug combinations
obs_targets <- obs_targets %>%
  mutate(perturb = case_when(
    drug == "Gefitinib_AZD4547_1to80" ~ paste(unique(obs_targets$perturb[obs_targets$drug == "Gefitinib"],
                                                     obs_targets$perturb[obs_targets$drug == "AZD4547"]), sep = ";"),
    drug == "LapatinibAZD4547" ~ paste(unique(obs_targets$perturb[obs_targets$drug == "Lapatinib"],
                                              obs_targets$perturb[obs_targets$drug == "AZD4547"]), sep = ";"),
    drug == "Selumetinib_MK2206_1to2" ~ paste(unique(obs_targets$perturb[obs_targets$drug == "Selumetinib"],
                                                     obs_targets$perturb[obs_targets$drug == "MK2206"]), sep = ";"),
    drug == "Selumetinib_MK2206_3to1" ~ paste(unique(obs_targets$perturb[obs_targets$drug == "Selumetinib"],
                                                     obs_targets$perturb[obs_targets$drug == "MK2206"]), sep = ";"),
    !drug %in% c("Gefitinib_AZD4547_1to80", "LapatinibAZD4547", "Selumetinib_MK2206_1to2", "Selumetinib_MK2206_3to1") ~ perturb
  )) %>%
  mutate(targets = case_when(
    drug == "Gefitinib_AZD4547_1to80" ~ paste(unique(obs_targets$targets[obs_targets$drug == "Gefitinib"],
                                                     obs_targets$targets[obs_targets$drug == "AZD4547"]), sep = ";"),
    drug == "LapatinibAZD4547" ~ paste(unique(obs_targets$targets[obs_targets$drug == "Lapatinib"],
                                              obs_targets$targets[obs_targets$drug == "AZD4547"]), sep = ";"),
    drug == "Selumetinib_MK2206_1to2" ~ paste(unique(obs_targets$targets[obs_targets$drug == "Selumetinib"],
                                                     obs_targets$targets[obs_targets$drug == "MK2206"]), sep = ";"),
    drug == "Selumetinib_MK2206_3to1" ~ paste(unique(obs_targets$targets[obs_targets$drug == "Selumetinib"],
                                                     obs_targets$targets[obs_targets$drug == "MK2206"]), sep = ";"),
    !drug %in% c("Gefitinib_AZD4547_1to80", "LapatinibAZD4547", "Selumetinib_MK2206_1to2", "Selumetinib_MK2206_3to1") ~ targets
  )) %>%
  rename("targets1" = perturb)

# Set perturb
if (perturb == "kinomebeads"){
  obs_targets <- obs_targets %>%
    rename("perturb" = targets)
} else {
  obs_targets <- obs_targets %>%
    rename("perturb" = targets1)
}

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

# filter out experiments where no activity was inferred for perturbed kinases
df_filtered <- map_dfr(1:nrow(df), function(i){
  tmp <- df[i,]
  targets <- target_df %>%
    filter(sample %in% tmp$experiment) %>%
    pull(perturb) %>%
    str_split(";") %>%
    unlist

  sum_act <- sum(tmp[,colnames(tmp) %in% targets], na.rm = T)

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

