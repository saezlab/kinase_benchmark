if(exists("snakemake")){
  curve_file <- snakemake@input$curve_file
  meta_file <- snakemake@input$meta_file
  mapped_file <- snakemake@input$mapped_file
  phospho_output <- snakemake@output$phospho
  meta_output <- snakemake@output$meta_output
}else{
  curve_file <- "data/decryptm/10_Kinase_Inhibitors/Phosphoproteome/curves_6KI.txt"
  meta_file <- "data/decryptm/10_Kinase_Inhibitors/Phosphoproteome/pipeline_6KI.toml"
  mapped_file <- "results/decryptm/protein_mapping/mapped_protein_6KI.csv"
  phospho_output <- "results/decryptm/phosphoproteome/intensities_6KI.csv"
  meta_output <- "results/decryptm/phosphoproteome/metadata_6KI.csv"
}

library(tidyverse)
library(configr)

## load files -------------------------
maxquant <- read.csv(curve_file, sep = "\t")
mapped_pps <- read.csv(mapped_file, sep = "\t")

# add protein position and fifteenmer to maxquant file
maxquant <- maxquant %>%
  left_join(mapped_pps %>%
              select(All.Phospho..STY..Probabilities, psite_location, peptide_start, fifteenmer),
            by = c("All.Phospho..STY..Probabilities"),
            relationship = "many-to-many") %>%
  mutate(pps_id = paste(Leading.proteins, Gene.names, psite_location, fifteenmer, sep = "|"))

## load and process meta data -------------------------
if (length(meta_file) != 1){
  meta_df <- map_dfr(meta_file, function(meta){
    metadata <- read.config(file = meta)
    drugs <- str_split(metadata$`Meta Data`$experiment_description, ", ") %>%
      unlist()

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(channels = as.character(metadata$TMT$channels),
                 drug = drug,
                 dose = metadata$TMT$doses,
                 dose_scale = metadata$`Meta Data`$dataset_name,
                 dataset = metadata$`Meta Data`$cells,
                 time = metadata$`Meta Data`$treatment_time,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 experiment_type = metadata$`Meta Data`$experiment_type) %>%
        mutate(drug = recode(drug,
                             "Rafametinib" = "Refametinib",
                             "Stausporin" = "Staursporin"))  %>%
        mutate(control = channels == metadata$Processing$control_channel) %>%
        mutate(time_min = str_remove(time, " min")) %>%
        mutate(sample = case_when(
          !control ~ paste0(experiment_type, "PTM_", cells, "_", drug, "_", time_min, "min_", channels, "_", dose),
          control ~ paste0(experiment_type, "PTM_", cells, "_", drug, "_", time_min, "min_", channels, "_", dose, "C")
        )) %>%
        mutate(Experiment = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time_min, "min_R1"))
    })
  })

} else {
  metadata <- read.config(file = meta_file)
  drugs <- str_split(metadata$`Meta Data`$experiment_description, ", ") %>%
    unlist()

  meta_df <- map_dfr(drugs, function(drug){
    data.frame(channels = as.character(metadata$TMT$channels),
               drug = drug,
               dose = metadata$TMT$doses,
               dose_scale = metadata$TMT$dose_scale,
               dataset = metadata$`Meta Data`$dataset_name,
               time = metadata$`Meta Data`$treatment_time,
               cells = metadata$`Meta Data`$cells,
               who = metadata$`Meta Data`$who,
               experiment_type = metadata$`Meta Data`$experiment_type) %>%
      mutate(drug = recode(drug,
                           "Rafametinib" = "Refametinib",
                           "Stausporin" = "Staursporin")) %>%
      mutate(control = channels == metadata$Processing$control_channel) %>%
      mutate(time_min = str_remove(time, " min")) %>%
      mutate(sample = case_when(
        !control ~ paste0(experiment_type, "PTM_", cells, "_", drug, "_", time_min, "min_", channels, "_", dose),
        control ~ paste0(experiment_type, "PTM_", cells, "_", drug, "_", time_min, "min_", channels, "_", dose, "C")
      )) %>%
      mutate(Experiment = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time_min, "min_R1"))
  })
}

head(meta_df)

## Format phosphorylation sites into matrix -------------------------
# select reporter intensities, phosphorylation sites and experiments
psp <- colnames(maxquant)[str_detect(colnames(maxquant), "Reporter.intensity.corrected")]
psp <- c(psp, c("pps_id", "Experiment"))

df <- maxquant %>%
  filter(!psite_location == "") %>%
  select(all_of(psp)) %>%
  mutate(tmp_id = paste(pps_id, Experiment, sep = ":")) %>% # create tmp ID to check for pps that were measured twice in each experiment
  select(-c(pps_id, Experiment)) %>%
  distinct() # remove duplicated rows

# remove rows with only zeros
msk <- rowSums(df[,1:11]) == 0
df_filtered <- df[!msk,]

# identify duplicated phosphorylation sites per experiment (e.g. due to missed cleavage)
dup_peptides <- df_filtered %>%
  filter(duplicated(tmp_id)) %>%
  pull(tmp_id) %>%
  unique()
length(dup_peptides)

# For each duplicated pps check correlation
# Select the one with higher intensities for pps with a high correlation
# remove the ones with a low correlation
fixed_dup <- map_dfr(dup_peptides, function(dup){
  df_dup <- df_filtered %>% filter(tmp_id == dup)
  pear_cor <- cor(df_dup[1,names(df_dup)[!names(df_dup) %in% c("pps_id", "tmp_id")]] %>% as.numeric(),
                  df_dup[2,names(df_dup)[!names(df_dup) %in% c("pps_id", "tmp_id")]] %>% as.numeric(),
                  method = "pearson")

  if (pear_cor >= 0.8){
    idx <- which(c(df_dup[1,1], df_dup[2,1]) == max(c(df_dup[1,1], df_dup[2,1])))[1]

    df_dup[idx,] %>% add_column(correlation = pear_cor)
  } else {
    df_dup[1,] %>% mutate(tmp_id = "remove") %>% add_column(correlation = pear_cor)
  }
}) %>%
  filter(tmp_id != "remove") %>%
  select(-correlation)

length(unique(fixed_dup$tmp_id))

# remove duplicates
df_rm_dup <- df_filtered %>%
  filter(!tmp_id %in% dup_peptides) %>%
  rbind(fixed_dup)

# final check for duplicates
df_rm_dup$tmp_id %>% duplicated() %>% any()

df_long <- df_rm_dup %>%
  pivot_longer(!tmp_id, names_to = "channel", values_to = "intensity") %>%
  mutate(pps_id = map_chr(str_split(tmp_id, ":"), 1)) %>%
  mutate(Experiment = map_chr(str_split(tmp_id, ":"), 2)) %>%
  mutate(channels = str_remove(channel, "Reporter.intensity.corrected."))

# add channel information from meta data
df_long <- df_long %>%
  left_join(meta_df %>%
              select(channels, dose, control, Experiment, sample) %>%
              distinct(),
            by = c("channels", "Experiment"), relationship = "many-to-many")
pps_m <- df_long %>%
  select(pps_id, sample, intensity) %>%
  pivot_wider(names_from = sample, values_from = intensity, values_fill = NA)

write_csv(pps_m, phospho_output)
write_csv(meta_df, meta_output)
