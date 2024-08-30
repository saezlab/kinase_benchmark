if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  PKN <- snakemake@input$file_PKN
  PKN_name <- snakemake@wildcards$PKN
  csv <- snakemake@output$rds
  log <- snakemake@output$log
  output_folder <- snakemake@params$output_folder
  minsize <- snakemake@params$minsize
}else{
  dataset <- "results/01_processed_data/hernandez/datasets/benchmark.gct"
  PKN <- "results/01_processed_data/hernandez/mapped_priors/ptm-sea/phosphositeplus_phosformer15.gmt"
  PKN_name <- "phosphositeplus_phosformer15"
  output_folder <- "results/02_activity_scores/hernandez/ptmsea"
  log <- "results/02_activity_scores/hernandez/ptmsea/log/phosphositeplus_phosformer15.log"
  if(!require("ssGSEA2")) remotes::install_github('smuellerd/ssGSEA2', upgrade='never')
  minsize <- 1
}
minsize <- as.double(minsize)

## Libraries ---------------------------
library(ssGSEA2)
library(tidyselect)
library(tidyverse)

## Run ssGSEA/PTM-SEA ---------------------------
net <- read.csv(PKN, header = F)

map_dfr(1:nrow(net), function(row_ids){
  targets <- net[row_ids,] %>% str_split("\t") %>% unlist()
  data.frame(kin = targets[2], length = length(targets) - 2)
}) %>% dplyr::filter(length >= 100000)

res <- run_ssGSEA2(dataset,
                  output.prefix = PKN_name,
                  gene.set.databases = PKN,
                  output.directory = output_folder,
                  sample.norm.type = "none",
                  weight = 0.75,
                  correl.type = "rank",
                  statistic = "area.under.RES",
                  output.score.type = "NES",
                  nperm = 100,
                  min.overlap = minsize,
                  max.overlap = 100000,
                  extended.output = TRUE,
                  global.fdr = FALSE,
                  log.file = log)

ptmsea <- read.delim(file=paste0(output_folder, "/", PKN_name, "-scores.gct"), skip=2)
ptmsea <- ptmsea %>%
  dplyr::select(colnames(ptmsea)[!str_detect(colnames(ptmsea), "Signature")]) %>%
  dplyr::select(-No.columns.scored)
colnames(ptmsea) <- str_remove(colnames(ptmsea), "^X")

ptmsea <- ptmsea %>%
  pivot_longer(!id, names_to = "condition", values_to = "score") %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::rename(source = id) %>%
  add_column(method = "ptmsea") 

write_csv(ptmsea, csv)


