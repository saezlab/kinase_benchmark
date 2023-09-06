if(exists("snakemake")){
  ptmsig_file <- snakemake@input$ptmsig
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ptmsig_file <- "data/prior/ptm.sig.db.all.uniprot.human.v1.9.0.gmt"
  output_file <- "results/hernandez/prior/ptmsigdb.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
}

## Libraries ---------------------------
library(biomaRt)
library(qusage)
library(tidyverse)

## Construct kinase-substrate interaction network ---------------------------
PTMsigDB <- qusage::read.gmt(ptmsig_file)
PTMsigDB.kinase <- PTMsigDB[str_detect(names(PTMsigDB), "KINASE")]

PTMsig_df <- map_dfr(names(PTMsigDB.kinase), function(kin){
  data.frame(source = kin,
             target_site = PTMsigDB.kinase[[kin]])
})

PTMsig_df <- PTMsig_df %>%
  mutate(source = str_remove(source, "KINASE-PSP_")) %>%
  mutate(target_site = str_remove(target_site, "-p;u")) %>%
  mutate(target = map_chr(str_split(map_chr(str_split(target_site, ";"), 1), "-"), 1)) %>%
  mutate(position = map_chr(str_split(target_site, ";"), 2))

## Change kinases to common gene names---------------------------
PTMsig_df$source <- map_chr(str_split(paste(PTMsig_df$source,
                                            PTMsig_df$source,
                                            sep = "/"),
                                      "/"),
                            2)

# Add uniprot id of kinases to identify autophosphorylation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res_kin <- getBM(attributes = c('uniprot_gn_id',
                                'external_gene_name'),
                 values = PTMsig_df$target,
                 mart = mart)  %>%
  dplyr::rename("target_gene_name" = external_gene_name) %>%
  dplyr::rename("target" = uniprot_gn_id)

PTMsig_df <- left_join(PTMsig_df, res_kin, by = "target", relationship = "many-to-many")

PTMsig_df <- PTMsig_df%>%
  mutate(target_site = paste(target_gene_name, position, sep = "_"))

## Prepare data ---------------------------
hernandez_df <- read_csv(dataset)
identifiers <- hernandez_df %>%
  dplyr::select(ID) %>%
  dplyr::mutate(gene_name = map_chr(str_split(ID, "\\|"), 1))  %>%
  dplyr::mutate(position = map_chr(str_split(ID, "\\|"), 2)) %>%
  dplyr::mutate(ensg = map_chr(str_split(ID, "\\|"), 3)) %>%
  dplyr::mutate(ensp = map_chr(str_split(ID, "\\|"), 4)) %>%
  dplyr::mutate(target_site = paste0(gene_name, "_", position)) %>%
  dplyr::rename("site" = ID)


## Merge with network ---------------------------
PTMsig_prior <- left_join(PTMsig_df,
                          identifiers, by = "target_site", relationship = "many-to-many")

PTMsig_prior_df <- PTMsig_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  mutate(target = case_when(
    str_detect(pattern = source, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = source, string = target_site) ~ target
  )) %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(source, target, mor)


# Remove duplicated edges (if present)
PTMsig_prior_df <- PTMsig_prior_df %>%
  distinct()

## Save processed PTMsigDB
write_tsv(PTMsig_prior_df, output_file)
