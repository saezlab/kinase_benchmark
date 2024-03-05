#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input
  output_file <- snakemake@output
}else{
  fifteenmer_file <- "results/hernandez/processed_data/fifteenmer.csv"
  output_file <- "results/hernandez/processed_data/pps_fifteenmer.csv"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
fifteenmer <- read_csv(fifteenmer_file)

fifteenmer_ids <- fifteenmer %>%
  mutate(pps = paste0(Protein, "_", Aminoacid, Position))


fifteenmer_ids$sequence_corrected <- map_chr(1:nrow(fifteenmer_ids), function(r_idx){
  df_i <- fifteenmer_ids[r_idx,]
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

merged_df <- fifteenmer_ids %>%
  dplyr::select(Protein, Position, Aminoacid, sequence_corrected) %>%
  filter(!is.na(sequence_corrected)) %>%
  dplyr::rename("Sequence" = sequence_corrected)

write_csv(merged_df, output_file)
