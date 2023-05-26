if(exists("snakemake")){
  networkin_file <- snakemake@input$networkin_file
  output_file <- snakemake@output$tsv
  file_datasets <- snakemake@input$file_dataset
  networkin_score <- snakemake@input$score
  GPS_file <- snakemake@input$GPS_file
  output_file_merge <- snakemake@output$tsv_merge
}else{
  networkin_file <- "data/prior/networkin_human_predictions_3.1.tsv"
  output_file <- "results/prior/networkin.tsv"
  GPS_file <- "results/prior/GPS.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
  networkin_score <- 5
  output_file_merge <- "results/prior/GPSnetworkin.tsv"
}
networkin_score <- as.numeric(networkin_score)

## Libraries ---------------------------
library(tidyverse)

## load NetworKIN ---------------------------
# load NetworKin prediction scores; convert kinase names to gene symbol and filter to kinases in KS table
nwkin <- read.table(networkin_file, sep = "\t",  stringsAsFactors = F, header= T, comment.char = "")

# rename kinases
nwkin$id <- toupper(nwkin$id)
nwkin <- nwkin %>%
  mutate(id = recode(id,
                     "ABL" = "ABL1",
                     "AURORAA" = "AURKA",
                     "AURORAA" = "AURKA",
                     "AURORAB" = "AURKB",
                     "AURORAC" = "AURKC",
                     "BRK" = "PTK6",
                     "CAMKIALPHA" = "CAMK1",
                     "CAMKIIALPHA" = "CAMK2A",
                     "CAMKIIBETA" = "CAMK2B",
                     "CAMKIIDELTA" = "CAMK2D",
                     "CAMKIIGAMMA" = "CAMK2G",
                     "CAMKIV" = "CAMK4",
                     "CAMKIDELTA" = "CAMK1D",
                     "CK1ALPHA" = "CSNK1A1",
                     "CK1DELTA" = "CSNK1D",
                     "CK1EPSILON" = "CSNK1E",
                     "CK1GAMMA1" = "CSNK1G1",
                     "CK1GAMMA2" = "CSNK1G2",
                     "CK1GAMMA3" = "CSNK1G3",
                     "CK2ALPHA" = "CSNK2A1",
                     "CK2A2" = "CSNK2A2",
                     "DNAPK" = "PRKDC",
                     "GSK3ALPHA" = "GSK3A",
                     "GSK3BETA" = "GSK3B",
                     "IKKALPHA" = "CHUK",
                     "IKKBETA" = "IKBKB",
                     "LKB1" = "STK11",
                     "MRCKA" = "CDC42BPA",
                     "MST2" = "STK3",
                     "MRCKA" = "CDC42BPA",
                     "PDGFRBETA" = "PDGFRB",
                     "PDGFRALPHA" = "PDGFRA",
                     "PDHK1" = "PDK1",
                     "PDHK2" = "PDK2",
                     "PDHK3" = "PDK3",
                     "PKBALPHA" = "AKT1",
                     "PKBBETA" = "AKT2",
                     "PKBGAMMA" = "AKT3",
                     "PKCALPHA" = "PRKCA",
                     "PKCEPSILON" = "PRKCE",
                     "PKCETA" = "PRKCH",
                     "PKCGAMMA" = "PRKCG",
                     "PKCIOTA" = "PRKCI",
                     "PKCTHETA" = "PRKCQ",
                     "PKG1CGKI" = "PRKG1",
                     "TRKA" = "NTRK1",
                     "TRKB" = "NTRK2",
                     "YES" = "YES1",
                     "PKAALPHA" = "PRKACA",
                     "PKABETA" = "PRKACB",
                     "PKAGAMMA" = "PRKACG"))

# Select targets above certain NetworKIN score
nwkin_filtered <- nwkin[nwkin$networkin_score >= networkin_score, ]

## Mapping ---------------------------
# use the 11mer combined with target gene name for NetworKin
# matching with gene and 15mer in CPTAC data)
nwkin_filtered$sequence <- toupper(nwkin_filtered$sequence)
nwkin_filtered$sequence <- gsub("\\-", "_", nwkin_filtered$sequence)

# Prepare Site in NetworKIN
nwkin_filtered$Site1 <- paste0(nwkin_filtered$substrate_name, "|", nwkin_filtered$sequence)

#map this id to Uniprot/Swissprot id used for GPS GS set; first establish ensembl > Swissprot key
mappings <- read.table("data/knowledge_based_isoform_selection_v1.4.txt", sep = "\t", stringsAsFactors = F, header = T)
#mappings <- mappings[!duplicated(paste0(mappings$gene, mappings$protein, mappings$SwissProt)), ]
mappings_sp <- mappings[!is.na(mappings$SwissProt), ]
mappings_sp <- mappings_sp[!duplicated(paste0(mappings_sp$protein, mappings_sp$SwissProt)), ]
mappings_sp$protein <- sub("\\..*", "", mappings_sp$protein)
rownames(mappings_sp) <- mappings_sp$protein

#use gene symbol instead of SwissProt id
mappings_hugo <- mappings[!duplicated(paste0(mappings$gene, mappings$gene_name_BCM_refined)), ]
rownames(mappings_hugo) <- mappings_hugo$gene


# Prepare phosphorylation sites in data sets
pps <- map_dfr(file_datasets, function(file){
  df <- read_tsv(file, col_types = cols())
  data.frame(site = df$site)
})

pps <- pps %>%
  distinct(.keep_all = TRUE)

# Change fifteenmer to elevenmer for mapping
pps$fifteenmer <- sub("^[^|]*\\|[^|]*\\|[^|]*\\|", "", pps$site)
pps$fifteenmer <- sub("\\|.*", "", pps$fifteenmer)
pps$gene <- mappings_hugo[sub("\\|.*", "", pps$site), "gene_name_BCM_refined"]
pps$Site1 <- paste0(pps$gene, "|", sub("\\D\\D$","", sub("^\\D\\D", "",pps$fifteenmer)))
site_check <- pps$fifteenmer[duplicated(pps$fifteenmer)]


# Map pps from datasets to NetworKIN
nwkin_df <- left_join(nwkin_filtered, pps, by = "Site1", relationship = "many-to-many") %>%
  mutate(target = case_when(
    is.na(site) ~ Site1,
    !is.na(site) ~ site
  )) %>%
  mutate(target = case_when(
    id == substrate_name ~ paste0(target, "|auto"), #mark autophosphorylation
    id != substrate_name ~ target
  )) %>%
  dplyr::rename("source" = id) %>%
  dplyr::select(source, target) %>%
  add_column(mor = 1) %>%
  distinct()


## Save prior ---------------------------
write_tsv(nwkin_df, output_file)

## Combine NetworKIN and GPS gold standard ---------------------------
GPS <- read_tsv(GPS_file, col_types = cols())

## only add pps from kinases already present in GPS
nwkin_df_kin <- nwkin_df %>%
  dplyr::filter(source %in% GPS$source)

merge_df <- rbind(GPS, nwkin_df_kin)
merge_df <- merge_df[!duplicated(merge_df),]

## Save merged ---------------------------
write_tsv(merge_df, output_file_merge)
