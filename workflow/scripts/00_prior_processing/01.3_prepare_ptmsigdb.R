if(exists("snakemake")){
  ptmsig_file <- snakemake@input$ptmsig
  output_file <- snakemake@output$tsv
}else{
  ptmsig_file <- "data/kinase_libraries/prior/ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
  output_file <- "results/00_prior/raw/ptmsigdb.tsv"
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
  mutate(source = str_remove(source, "KINASE-iKiP_")) %>%
  mutate(target_site = str_remove(target_site, "-p;u")) %>%
  mutate(target = map_chr(str_split(map_chr(str_split(target_site, ";"), 1), "-"), 1)) %>%
  mutate(position = map_chr(str_split(target_site, ";"), 2))

## Change kinases to common gene names---------------------------
PTMsig_df$source <- map_chr(str_split(paste(PTMsig_df$source,
                                            PTMsig_df$source,
                                            sep = "/"),
                                      "/"),
                            2)

# Translate target uniprot id to gene name
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
PTMsig_prior_df <- PTMsig_df %>%
  dplyr::select(source, target_site) %>%
  dplyr::rename("target" = target_site) %>%
  mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
  mutate(position = map_chr(str_split(target, "_"), 2)) %>%
  mutate(mor = 1) %>%
  dplyr::select(source, target, target_protein, position, mor) %>%
  distinct() %>%
  mutate(sequence = NA)

## Save processed PTMsigDB
write_tsv(PTMsig_prior_df, output_file)
