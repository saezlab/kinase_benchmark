require(Hmisc)

calculate_Kinase_Activity <- function(mat, network, col="source", col2="target", min_sites=5, samps=colnames(mat)){
  #samps <- samps[!samps %in% c("fifteenmer", "prot", "sp", "site", "gene")]
  mat <- data.matrix(mat[ , samps])
  kins <- unique(unlist(network[ , col]))
  network <- network[!duplicated(paste0(network[ , col], network[, col2])), ]
  mean_scores <- matrix(NA, nrow = length(kins), ncol = ncol(mat), dimnames = list(kins, colnames(mat)))
  med_scores <- mean_scores
  uq_scores <- mean_scores
  pc1_scores <- mean_scores
  wilcox_scores <- mean_scores
  ks_scores <- mean_scores
  num_targets <- mean_scores

  for(i in 1:length(kins)){
    #print(i)
    targets <- intersect(rownames(mat), unlist(network[network[ , col]==kins[i], col2]))
    nontargets <- setdiff(rownames(mat), targets)
    if(length(targets) >= min_sites){
      mean_scores[kins[i], ] <- colMeans(mat[targets, ], na.rm = T)
      med_scores[kins[i], ] <- as.numeric(apply(mat[targets, ], 2, median, na.rm=T))
      uq_scores[kins[i], ] <- as.numeric(apply(mat[targets, ], 2, quantile, probs=0.75, na.rm=T))
      phos_imp <- mat[targets, ]
      for(j in 1:ncol(mat)){
        phos_imp[is.na(mat[targets, j]), j] <- mean_scores[kins[i], j]
      }
      phos_pca <- prcomp(t(phos_imp[ , colSums(is.na(phos_imp))==0]), scale=F)
      pc1_scores[kins[i], colSums(is.na(phos_imp))==0] <- t(phos_pca$x[ , 1])

      #also get number of targets in a given sample and set scores for methods that calculate scores without a minimum target requirement to NA if the minimum threshold is not met
      for(j in 1:ncol(mat)){
        target_num <- sum(!is.na(mat[targets, j]))
        num_targets[kins[i], j] <- target_num
        if(target_num >= min_sites){
          KS_p.val <- ks.test(as.numeric(mat[targets, j]), as.numeric(mat[nontargets, j]))$p.value
          if(KS_p.val == 0){
            KS_p.val <- 2.2e-16
          }
          ks_scores[kins[i], j] <- -log10(KS_p.val)
          wilcox_p.val <- wilcox.test(as.numeric(mat[targets, j]), as.numeric(mat[nontargets, j]))$p.value
          if(wilcox_p.val == 0){
            wilcox_p.val <- 2.2e-16
          }
          wilcox_scores[kins[i], j] <- -log10(wilcox_p.val)
          if(median(as.numeric(mat[targets, j]), na.rm = T) < median(as.numeric(mat[nontargets, j]), na.rm = T)){
            ks_scores[kins[i], j] <- -ks_scores[kins[i], j]
            wilcox_scores[kins[i], j] <- -wilcox_scores[kins[i], j]
          }
        } else {
          mean_scores[kins[i], j] <- NA
          pc1_scores[kins[i], j] <- NA
          med_scores[kins[i], j] <- NA
          uq_scores[kins[i], j] <- NA
        }
      }

      #for pca method, direction doesn't always track with direction of activity; compare PCA scores with mean scores to determine if direction should be adjusted
      if(sum(!is.na(mean_scores[kins[i], ])) > 4){
        if(rcorr(as.numeric(pc1_scores[kins[i], ]), as.numeric(mean_scores[kins[i], ]))$r[1, 2] < 0){
          pc1_scores[kins[i], ] <- -pc1_scores[kins[i], ]
        } else if(sum(!is.na(mean_scores[kins[i], ])) > 0){
          pos_corr <- sum(sign(mean_scores[kins[i], ])==sign(pc1_scores[kins[i], ]), na.rm = T)
          neg_corr <- sum(sign(mean_scores[kins[i], ])!=sign(pc1_scores[kins[i], ]), na.rm = T)
          if(neg_corr > pos_corr){
            pc1_scores[kins[i], ] <- -pc1_scores[kins[i], ]
          }
        }
      }

    } else {
      mean_scores <- mean_scores[rownames(mean_scores) != kins[i], ]
      med_scores <- med_scores[rownames(med_scores) != kins[i], ]
      uq_scores <- uq_scores[rownames(uq_scores) != kins[i], ]
      pc1_scores <- pc1_scores[rownames(pc1_scores) != kins[i], ]
      ks_scores <- ks_scores[rownames(ks_scores) != kins[i], ]
      wilcox_scores <- wilcox_scores[rownames(wilcox_scores) != kins[i], ]
      num_targets <- num_targets[rownames(num_targets) != kins[i], ]
    }
  }

  num_targets[is.na(num_targets)] <- 0
  all_scores <- list(mean_scores, med_scores, uq_scores, pc1_scores, ks_scores, wilcox_scores, num_targets)
  names(all_scores) <- c("mean","median","UQ","PC1","KS","Wilcox","number_of_targets")
  return(all_scores)
}
