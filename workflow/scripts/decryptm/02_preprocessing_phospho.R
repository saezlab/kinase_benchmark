if(exists("snakemake")){
  folder <- snakemake@input$input_folder
  mapped_file <- snakemake@input$mapped_file
  phospho_output <- snakemake@output$phospho
  meta_output <- snakemake@output$meta_output
}else{
  folder <- "data/decryptm/3_EGFR_Inhibitors"
  mapped_file <- "results/decryptm/protein_mapping/mapped_protein_3_EGFR_Inhibitors.csv"
  phospho_output <- "results/decryptm/phosphoproteome/EC50_3_EGFR_Inhibitors.csv"
  meta_output <- "results/decryptm/phosphoproteome/metadata_3_EGFR_Inhibitors.csv"
}

library(tidyverse)
library(configr)

## load files -------------------------
maxquant_files <- list.files(folder, pattern = "curve", recursive = T, full.names = T)
maxquant_files <- maxquant_files[str_detect(maxquant_files, "Phospho")]

meta_file <- list.files(folder, pattern = "toml", recursive = T, full.names = T)
meta_file <- meta_file[str_detect(meta_file, "Phospho")]

maxquant <- map_dfr(maxquant_files, function(curve_file) {
  df <- read.csv(curve_file, sep = "\t")
  df <- df[!colnames(df) %in% c("Phosphoproteome", "Fullproteome")]
  if (folder == "data/decryptm/Dasatinib_Triplicates"){
    file_id <- str_extract(curve_file, "(?<=Phosphoproteome_).*?(?=\\/curves.txt)")
  } else  if (folder == "data/decryptm/BCRABL_Inhibitors"){
    file_id <- "BCRABL"
  } else {
    file_id <- str_extract(curve_file, "(?<=curves_).*?(?=\\.txt)")
  }
  if(is.na(file_id)){
    file_id <- ""
  }
  df %>% mutate(sample = paste(Experiment, file_id, sep = "_"))
})
mapped_pps <- read.csv(mapped_file, sep = "\t")

# add protein position and fifteenmer to maxquant file
maxquant <- maxquant %>%
  left_join(mapped_pps %>%
              select(All.Phospho..STY..Probabilities, psite_location, peptide_start, fifteenmer),
            by = c("All.Phospho..STY..Probabilities"),
            relationship = "many-to-many") %>%
  mutate(pps_id = paste(Leading.proteins, Gene.names, psite_location, fifteenmer, sep = "|"))

## load and process meta data -------------------------
# The processing has to be done seperately for the different experiments as they differ in their nomenclature
if(folder == "data/decryptm/BCRABL_Inhibitors"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- "BCRABL"
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- c("Imatinib", "Dasatinib")

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = "R1",
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", ""))
    })
  })
} else if (folder == "data/decryptm/Dasatinib_Triplicates"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- str_split(metadata$`Meta Data`$experiment_description, ", ") %>%
      unlist()

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = c("Dasatinib"),
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1", "R2", "R3"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = "") %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", ""))
    })
  })
  meta_df$dataset_id <- rep(c("MS2", "MS3"), each = 3)
  meta_df <- meta_df %>%
    mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
    mutate(sample = str_replace(sample, " ", ""))
} else if (folder == "data/decryptm/Combination_Treatments_Selumetinib_MK2206"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- dataset_id

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1", "R2"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", "")) %>%
        mutate(sample = recode(sample,
                               "ddPTM_A459_Selumetinib_MK2206_1to2_30min_R1_Selumetinib_MK2206_1to2" = "ddPTM_A459_SelumetinibMK2206-1to2_30min_R1_Selumetinib_MK2206_12",
                               "ddPTM_A459_Selumetinib_MK2206_3to1_30min_R1_Selumetinib_MK2206_3to1" = "ddPTM_A459_SelumetinibMK2206-3to1_30min_R1_Selumetinib_MK2206_31"))
      })
  })
} else if (folder == "data/decryptm/Combination_Treatments_AZD4547_Lapatinib"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- c("Lapatinib", "AZD4547", "LapatinibAZD4547")

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1", "R2"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", ""))
    })
  })
} else if (folder == "data/decryptm/Combination_Treatments_AZD4547_Gefitinib"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- dataset_id

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", "")) %>%
        mutate(sample = recode(sample,
                               "ddPTM_PC-9_Gefitinib_AZD4547_1to80_30min_R1_Gefitinib_AZD4547_1to80" = "ddPTM_PC-9_GeftinibAZD4547-1to80_30min_R1_Gefitinib_AZD4547"))
    })
  })
} else if (folder == "data/decryptm/HER2_Inhibitors"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- str_split(dataset_id, "_")[[1]][2]

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = "2h",
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", "")) %>%
        mutate(cells = recode(cells,
                              "BTB-474" = "BT-474")) %>%
        mutate(sample = recode(sample,
                               "ddPTM_BTB-474_Trastuzumab_2h_R1_BT474_Trastuzumab" = "ddPTM_BT-474_Trastuzumab_2h_R1_BT474_Trastuzumab"))
      })
  })
} else if (folder == "data/decryptm/10_Kinase_Inhibitors"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- str_split(metadata$`Meta Data`$experiment_description, ", ") %>%
      unlist()

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(drug = recode(drug,
                             "Rafametinib" = "Refametinib",
                             "Stausporin" = "Staursporin")) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", ""))
    })
  })
} else if (folder == "data/decryptm/3_EGFR_Inhibitors"){
  meta_df <- map_dfr(meta_file, function(meta){
    dataset_id <- str_extract(meta, "(?<=pipeline_).*?(?=\\.toml)")
    if(is.na(dataset_id)){
      dataset_id <- ""
    }
    metadata <- read.config(file = meta)
    drugs <- c("Afatinib", "Gefitinib", "Dasatinib")

    meta_df <- map_dfr(drugs, function(drug){
      data.frame(drug = drug,
                 dose = paste(metadata$TMT$doses, collapse = ";"),
                 channels = paste(metadata$TMT$channels, collapse = ";"),
                 control_channel = metadata$Processing$control_channel,
                 dose_scale = metadata$TMT$dose_scale,
                 dose_unit = metadata$TMT$dose_label,
                 dataset = metadata$`Meta Data`$dataset_name,
                 time = metadata$`Meta Data`$treatment_time,
                 quantification = metadata$TMT$quantification,
                 cells = metadata$`Meta Data`$cells,
                 who = metadata$`Meta Data`$who,
                 replicate = c("R1", "R2", "R3", "R4"),
                 experiment_type = metadata$`Meta Data`$experiment_type,
                 dataset_id = dataset_id) %>%
        mutate(drug = recode(drug,
                             "Rafametinib" = "Refametinib",
                             "Stausporin" = "Staursporin")) %>%
        mutate(sample = paste0(experiment_type, "PTM_", cells, "_", drug, "_", time, "_", replicate, "_", dataset_id)) %>%
        mutate(sample = str_replace(sample, " ", ""))
    })
  })
}

maxquant$Experiment %>% unique()

## Format phosphorylation sites into matrix -------------------------
# select EC50 values, phosphorylation sites and experiments
psp <- c("pps_id", "sample", "R2", "pEC50", "EC50", "Curve.effect.size")

df <- maxquant %>%
  filter(!psite_location == "") %>%
  select(all_of(psp)) %>%
  mutate(tmp_id = paste(pps_id, sample, sep = ":")) %>% # create tmp ID to check for pps that were measured twice in each experiment
  #filter(!is.na(Log.EC50)) %>%
  select(-c(pps_id, sample)) %>%
  distinct() # remove duplicated rows

# remove rows with only zeros
msk <- rowSums(df[c("R2", "pEC50", "EC50", "Curve.effect.size")], na.rm = T) == 0
df <- df[!msk,]

ggplot(df, aes(x = R2)) +
  geom_histogram(bins = 50)

# remove pps with R2 of 0
df <- df %>%
  filter(!R2 == 0)

# identify duplicated phosphorylation sites per experiment (e.g. due to missed cleavage)
dup_peptides <- df %>%
  filter(duplicated(tmp_id)) %>%
  pull(tmp_id) %>%
  unique()
length(dup_peptides)

# For each duplicated pps check correlation of intensities
int_rep <- colnames(maxquant)[str_detect(colnames(maxquant), "Reporter.intensity.corrected")]
int <- c(int_rep, c("pps_id", "sample"))

df_int <- maxquant %>%
  filter(!psite_location == "") %>%
  select(all_of(int)) %>%
  mutate(tmp_id = paste(pps_id, sample, sep = ":")) %>% # create tmp ID to check for pps that were measured twice in each experiment
  select(-c(pps_id, sample)) %>%
  distinct() # remove duplicated rows

# remove rows with only zeros
msk <- rowSums(df_int[int_rep], na.rm = T) == 0
df_int <- df_int[!msk,]

# Select the one with lower R2 for pps with a high correlation
# remove the ones with a low correlation
fixed_dup <- map_dfr(dup_peptides, function(dup){
  df_dup <- df_int %>% filter(tmp_id == dup)
  df_ec50 <- df %>% filter(tmp_id == dup)
  pear_cor <- cor(df_dup[1,names(df_dup)[!names(df_dup) %in% c("pps_id", "tmp_id")]] %>% as.numeric(),
                  df_dup[2,names(df_dup)[!names(df_dup) %in% c("pps_id", "tmp_id")]] %>% as.numeric(),
                  method = "pearson", use = "complete.obs")

  if (pear_cor >= 0.8){
    # Select phosphorylation site with higher intensities as these can be considered to be more reliable
    idx <- which(c(rowSums(df_dup[!colnames(df_dup) == "tmp_id"])) == max(c(rowSums(df_dup[!colnames(df_dup) == "tmp_id"]))))[1]
    #idx <- which(c(df_ec50$R2) == max(c(df_ec50$R2), na.rm = T))[1]

    df_ec50[idx,] %>% add_column(correlation = pear_cor)
  } else {
    df_ec50[1,] %>% mutate(tmp_id = "remove") %>% add_column(correlation = pear_cor)
  }
})

fixed_dup <- fixed_dup %>%
  filter(tmp_id != "remove") %>%
  select(-correlation)

length(unique(fixed_dup$tmp_id))

# remove duplicates
df_rm_dup <- df %>%
  filter(!tmp_id %in% dup_peptides) %>%
  rbind(fixed_dup)

# final check for duplicates
df_rm_dup$tmp_id %>% duplicated() %>% any()

pps_long <- df_rm_dup %>%
  mutate(sample = map_chr(str_split(tmp_id, ":"), 2)) %>%
  mutate(pps_id = map_chr(str_split(tmp_id, ":"), 1)) %>%
  select(-tmp_id)

write_csv(pps_long, phospho_output)
write_csv(meta_df, meta_output)
