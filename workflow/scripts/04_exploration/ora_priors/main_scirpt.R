library(piano)
library(GSEABase)
library(tidyverse)
library(reshape2)
library(pheatmap)

source("workflow/scripts/ora_priors/support_functions.R")

pathways_df <- data.frame(import_gmt("data/ora/c2.cp.v2023.1.Hs.symbols.gmt"))
phosphositeplus <- as.data.frame(read_delim("results/prior/phosphositeplus.tsv",
                              delim = "\t", escape_double = FALSE,
                              trim_ws = TRUE))

background <- unique(phosphositeplus$target_protein)
pathways_df <- pathways_df[which(pathways_df$gene %in% background),]

ora_res_list <- list()
for(kinase in unique(phosphositeplus$source))
{

  sig_set <- unique(phosphositeplus[which(phosphositeplus$source == kinase),"target_protein"])
  ora_res <- runGSAhyper(sig_set, universe = background, gsc = loadGSC(pathways_df), adjMethod = "fdr")

  res_table <-  as.data.frame(ora_res$resTab)
  res_table$kinase <- kinase

  ora_res_list[[kinase]] <- res_table
}

ora_res_df <- as.data.frame(do.call(rbind, ora_res_list))
ora_res_df$n_targets <- ora_res_df$`Significant (in gene set)` + ora_res_df$`Significant (not in gene set)`

ora_res_df_reduced <- ora_res_df[which(ora_res_df$n_targets > 20),]
ora_res_df_reduced$pathway <- gsub(".+[.]","",row.names(ora_res_df_reduced))

ora_res_df_reduced_wide <- dcast(data = ora_res_df_reduced_wide, formula = pathway~kinase, value.var = "p-value")
pheatmap(cor(ora_res_df_reduced_wide[,-1]), cutree_rows = 6, cutree_cols = 6)

regulation_power <- colSums(ora_res_df_reduced_wide[,-1] <= 0.05)
names(regulation_power) <- names(ora_res_df_reduced_wide[,-1])
regulation_power <- as.data.frame(regulation_power) %>%
  rownames_to_column("targets")

rank <- read_csv("results/manuscript_figures/figure_3/kinase_GSknown.csv")

df <- full_join(rank, regulation_power, by = "targets")
df %>%
  ggplot(aes(x = mean_rank, y = regulation_power)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.6, fullrange = TRUE, color = "steelblue3") +
  theme_bw()

lm(mean_rank ~ 1 + regulation_power, data = df)

plot(dens)

write_csv(ora_res_df, file = "results/ora_targets/ora_res_df.csv")
write_csv(ora_res_df_reduced_wide, file = "results/ora_targets/ora_res_df_reduced_wide.csv")
