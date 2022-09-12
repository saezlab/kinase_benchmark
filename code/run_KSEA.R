library("KSEAapp")

KSEA_input_files <- list.files(path = "output/CPTAC/KSEA", pattern = "input", full.names = T)

map(KSEA_input_files, function(file_path){
  KSEA_input <- read.csv(file_path)
  KS_network <- read.csv("data/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")

  KSEA_scores <- KSEA.Scores(KSData = KS_network,
                             PX = KSEA_input,
                             NetworKIN = F)

  write.csv(KSEA_scores, paste0(str_replace(file_path, "input", "output")), row.names = F)
})

