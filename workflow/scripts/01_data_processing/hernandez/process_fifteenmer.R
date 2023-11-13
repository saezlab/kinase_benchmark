#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input
  output_file <- snakemake@output
}else{
  input_file <- "results/hernandez/processed_data/benchmark_data.csv"
  fifteenmer_file <- "results/hernandez/processed_data/fifteenmer.csv"
  output_file <- "results/hernandez/processed_data/pps_fifteenmer.csv"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
data <- read_csv(input_file)
pps_ids <- data %>%
  select(ID) %>%
  mutate(protein = map_chr(str_split(ID, "\\|"), 1)) %>%
  mutate(aa = map_chr(str_split(ID, "\\|"), 2)) %>%
  mutate(pps = paste(protein, aa, sep = "_"))

fifteenmer <- read_csv(fifteenmer_file)

fifteenmer_ids <- fifteenmer %>%
  mutate(pps = paste0(Protein, "_", Aminoacid, Position))

merged <- left_join(pps_ids, fifteenmer_ids, by = "pps", relationship = "many-to-many")
merged <- merged %>%
  distinct(pps, .keep_all = TRUE)

merged$sequence_corrected <- map_chr(1:nrow(merged), function(r_idx){
  df_i <- merged[r_idx,]
  if(is.na(df_i$Sequence)) {
    new_sequence <- NA
  } else  {
    upstream <- map_chr(str_split(df_i$Sequence, "\\(ph\\)"), 1)
    downstream <- map_chr(str_split(df_i$Sequence, "\\(ph\\)"), 2)
    length_up <- as.numeric(nchar(upstream))
    length_down <- as.numeric(nchar(downstream))

    upstream_new <- paste0(rep("_", times = 8-length_up), collapse = "")
    downstream_new <- paste0(rep("_", times = 7-length_down), collapse = "")

    new_sequence <- paste0(upstream_new, upstream, downstream, downstream_new)
  }
})

nrow(merged %>%
       filter(!is.na(sequence_corrected)))/nrow(pps_ids)

merged_df <- merged %>%
  select(ID, pps, sequence_corrected) %>%
  filter(!is.na(sequence_corrected)) %>%
  rename("fifteenmer" = sequence_corrected)

write_csv(merged_df, output_file)
