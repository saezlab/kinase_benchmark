outlieR_plus <- function(input_df, reference_samples, testing_samples, z_thresh="default", z_method=NA, type="outlieR", scaling_factor=1.48, quant=0.95, genesAsRownames=T, gene_column=NULL, ...){
  OP_output <- list()
  if(!type %in% c("outlieR","COPA","OS","ORA","Tstat")){
    stop("type must be set to valid method: outlieR(default), COPA, OS, ORA, Tstat")
  }
  reference_samples <- intersect(reference_samples, colnames(input_df))
  testing_samples <- intersect(testing_samples, colnames(input_df))
  all_samples <- union(reference_samples, testing_samples)
  if(genesAsRownames==F){
    rownames(input_df) <- input_df[ , gene_column]
  }
  input_df <- input_df[ , all_samples]
  if(is.na(z_method)){
    if(type!="Tstat"){
      z_method <- "medMAD"
    } else {
      z_method <- "normal"
    }
  }
  if(z_method=="normal"){
    scaling_factor <- 1
  } 
  if(z_thresh=="default"){
    if(type %in% c("outlieR","Tstat","COPA")){
      z_thresh <- 2.3265
    } else {
    z_thresh <- 2.3265 * 1.48
    }
  }
  if(type=="outlieR"){
    # first calculate Z-scores for Z-score matrix
    OP_output[["Zscores"]] <- outlieR(input_df = input_df, ref_samples = reference_samples, test_samples = all_samples, z_method = z_method, z_thresh = z_thresh, scaling_factor=scaling_factor, ...)
  } else if(type=="Tstat"){
    OP_output[["Zscores"]] <- outlieR(input_df = input_df, ref_samples = all_samples, test_samples = all_samples, z_method = z_method, z_thresh = z_thresh, scaling_factor = scaling_factor, ...)
  } else if(type=="COPA"){
    OP_output[["Zscores"]] <- outlieR(input_df = input_df, ref_samples = all_samples, test_samples = all_samples, z_method = z_method, z_thresh = z_thresh, ...)
  } else if(type=="OS"){
    OP_output[["Zscores"]] <- outlieR(input_df = input_df, ref_samples = all_samples, scaling_sample_set1 = reference_samples, test_samples = all_samples, z_method = z_method, z_thresh = z_thresh, ...)
  } else if(type=="ORA"){
    OP_output[["Zscores"]] <- outlieR(input_df = input_df, ref_samples = all_samples, test_samples = all_samples, z_method = z_method, z_thresh = z_thresh, ...)
  } 
  
  # now identify outliers and generate outlier matrix (1 if outlier, 0 if not)  
  OP_output[["outliers"]] <- as.data.frame(matrix(NA, nrow = nrow(input_df), ncol = length(all_samples), dimnames = list(rownames(input_df), all_samples)))
  if(type=="outlieR"|type=="Tstat"){
    OP_output[["outliers"]] <- OP_output[["Zscores"]][ , all_samples] > z_thresh
    class(OP_output[["outliers"]]) <- "numeric"
    OP_output[["outliers"]] <- OP_output[["outliers"]] - (OP_output[["Zscores"]][ , all_samples] < -z_thresh)
  } else if(type=="COPA"){
    for(i in 1:nrow(OP_output[["outliers"]])){
      thresh1 <- quantile(OP_output[["Zscores"]][i, testing_samples], quant, na.rm = T)
      thresh2 <- quantile(OP_output[["Zscores"]][i, testing_samples], 1-quant, na.rm = T)
      OP_output[["outliers"]][i, ] <- as.numeric(OP_output[["Zscores"]][i, all_samples] > thresh1)
      OP_output[["outliers"]][i, ] <- OP_output[["outliers"]][i, ] - as.numeric(OP_output[["Zscores"]][i, all_samples] < thresh2)
    }
  } else if(type=="OS"){
    a <- 1
    for(i in 1:nrow(OP_output[["outliers"]])){
      Q1 <- as.numeric(quantile(OP_output[["Zscores"]][i, all_samples], 0.25, na.rm = T))
      Q3 <- as.numeric(quantile(OP_output[["Zscores"]][i, all_samples], 0.75, na.rm = T))
      iqr <- Q3 - Q1
      thresh1 <- (a * iqr) + Q3
      thresh2 <- Q1 - (a * iqr)
      OP_output[["outliers"]][i, ] <- as.numeric(OP_output[["Zscores"]][i, all_samples] > thresh1)
      OP_output[["outliers"]][i, ] <- OP_output[["outliers"]][i, ] - as.numeric(OP_output[["Zscores"]][i, all_samples] < thresh2)
    }
  } else if(type=="ORA"){
    a <- 1.5
    for(i in 1:nrow(OP_output[["outliers"]])){
      Q1 <- as.numeric(quantile(input_df[rownames(OP_output[["outliers"]])[i], all_samples], 0.25, na.rm = T))
      Q3 <- as.numeric(quantile(input_df[rownames(OP_output[["outliers"]])[i], all_samples], 0.75, na.rm = T))
      iqr <- Q3 - Q1
      thresh1 <- (a * iqr) + Q3
      thresh2 <- Q1 - (a * iqr)
      OP_output[["outliers"]][i, ] <- as.numeric(input_df[rownames(OP_output[["outliers"]])[i], all_samples] > thresh1)
      OP_output[["outliers"]][i, ] <- OP_output[["outliers"]][i, ] - as.numeric(input_df[rownames(OP_output[["outliers"]])[i], all_samples] < thresh2)
    }
  } 
  
  OP_output[["DEscores"]] <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(OP_output[["outliers"]]), dimnames = list(rownames(OP_output[["outliers"]]), c("PositiveOutlierScores","NegativeOutlierScores"))))
  if(type=="outlieR"){
    OP_output[["DEscores"]]$PositiveOutlierScores <- rowSums(OP_output[["outliers"]][ , testing_samples]==1, na.rm = T)/rowSums(!is.na(OP_output[["outliers"]][ , testing_samples])) - rowSums(OP_output[["outliers"]][ , reference_samples]==1, na.rm = T)/rowSums(!is.na(OP_output[["outliers"]][ , reference_samples]))
    OP_output[["DEscores"]]$NegativeOutlierScores <- rowSums(OP_output[["outliers"]][ , testing_samples]==(-1), na.rm = T)/rowSums(!is.na(OP_output[["outliers"]][ , testing_samples])) - rowSums(OP_output[["outliers"]][ , reference_samples]==(-1), na.rm = T)/rowSums(!is.na(OP_output[["outliers"]][ , reference_samples]))
  } else if(type=="OS"){
    Zout <- OP_output[["Zscores"]][ , colnames(OP_output[["outliers"]])] * as.numeric(OP_output[["outliers"]]==1)
    OP_output[["DEscores"]]$PositiveOutlierScores <- rowSums(Zout[ , testing_samples], na.rm = T)/rowSums(!is.na(OP_output[["outliers"]][ , testing_samples]))
    Zout <- OP_output[["Zscores"]][ , colnames(OP_output[["outliers"]])] * as.numeric(OP_output[["outliers"]]==(-1))
    OP_output[["DEscores"]]$NegativeOutlierScores <- rowSums(Zout[ , testing_samples], na.rm = T)/rowSums(!is.na(OP_output[["outliers"]][ , testing_samples]))
  } else if(type=="Tstat"){
    if(z_method=="normal"){
      Xm <- rowMeans(input_df[, reference_samples], na.rm = T)
      Ym <- rowMeans(input_df[ , testing_samples], na.rm = T)
      Xv <- apply(input_df[ , reference_samples], 1, function(x) sum((x-mean(x, na.rm=T))^2))
      Yv <- apply(input_df[ , testing_samples], 1, function(x) sum((x-mean(x, na.rm=T))^2))
      Ns <- rowSums(!is.na(input_df)) - 2
      OP_output[["DEscores"]]$PositiveOutlierScores <- (Ym - Xm)/sqrt((Xv + Yv)/(Ns + 0.0000000001) + 0.0000000001)
    } else {
      Xm <- apply(input_df[ , reference_samples], 1, median, na.rm = T)
      Ym <- apply(input_df[ , testing_samples], 1, median, na.rm = T)
      Xv <- apply(input_df[ , reference_samples], 1, function(x) abs(x-median(x, na.rm=T)))
      Yv <- apply(input_df[ , testing_samples], 1, function(x) abs(x-median(x, na.rm=T)))
      mads <- apply(cbind(t(Xv), t(Yv)), 1, median, na.rm=T)
      OP_output[["DEscores"]]$PositiveOutlierScores <- (Ym - Xm)/(mads + 0.0000000001)
    }
    OP_output[["DEscores"]]$NegativeOutlierScores <- OP_output[["DEscores"]]$PositiveOutlierScores
  } else if(type=="ORA"){
    nO <- cbind(rowSums(OP_output[["outliers"]][ , testing_samples]==1), rowSums(OP_output[["outliers"]]==1))
    nO[nO[ , 1]==0, 1] <- 1
    OP_output[["DEscores"]]$PositiveOutlierScores <- apply(nO, 1, function(x) phyper(x[1] - 1, length(testing_samples), length(all_samples) - length(testing_samples), x[2], lower.tail = F))
    OP_output[["DEscores"]]$PositiveOutlierScores <- -log10(OP_output[["DEscores"]]$PositiveOutlierScores)
    OP_output[["DEscores"]]$PositiveOutlierScores[nO[ , 2]==0] <- 0
    nO <- cbind(rowSums(OP_output[["outliers"]][ , testing_samples] == -1), rowSums(OP_output[["outliers"]] == -1))
    nO[nO[ , 1]==0, 1] <- 1
    OP_output[["DEscores"]]$NegativeOutlierScores <- apply(nO, 1, function(x) phyper(x[1] - 1, length(testing_samples), length(all_samples) - length(testing_samples), x[2], lower.tail = F))
    OP_output[["DEscores"]]$NegativeOutlierScores <- -log10(OP_output[["DEscores"]]$NegativeOutlierScores)
    OP_output[["DEscores"]]$NegativeOutlierScores[nO[ , 2]==0] <- 0
  } else if(type=="COPA"){
    OP_output[["DEscores"]]$PositiveOutlierScores <- apply(OP_output[["Zscores"]][ , testing_samples], 1, quantile, quant)
    OP_output[["DEscores"]]$NegativeOutlierScores <- apply(OP_output[["Zscores"]][ , testing_samples], 1, quantile, 1-quant)
  }
  #direction indicates whether or not the DE score for a given method is positive or negative for negative outliers
  if(type=="outlieR"|type=="ORA"){
    direction <- "pos"
  } else {
    direction <- "neg"
  }
  OP_output[["DEscores"]]$best <- OP_output[["DEscores"]]$PositiveOutlierScores
  #if direction is negative, set the best outlier score to the negative outlier score when the magnitude of this score is greater than the positive outlier score
  if(direction=="neg"){
    OP_output[["DEscores"]]$best[!is.na( OP_output[["DEscores"]]$best) & (-OP_output[["DEscores"]]$NegativeOutlierScores > OP_output[["DEscores"]]$PositiveOutlierScores)] <- OP_output[["DEscores"]]$NegativeOutlierScores[!is.na( OP_output[["DEscores"]]$best) & ((-OP_output[["DEscores"]]$NegativeOutlierScores) > OP_output[["DEscores"]]$PositiveOutlierScores)]
  } else if(direction=="pos"){
  #if direction is positive, set the best outlier score to -the negative outlier score when this score is greater than the positive outlier score
    OP_output[["DEscores"]]$best[!is.na( OP_output[["DEscores"]]$best) & (OP_output[["DEscores"]]$NegativeOutlierScores > OP_output[["DEscores"]]$PositiveOutlierScores)] <- -OP_output[["DEscores"]]$NegativeOutlierScores[!is.na( OP_output[["DEscores"]]$best) & ((OP_output[["DEscores"]]$NegativeOutlierScores) > OP_output[["DEscores"]]$PositiveOutlierScores)]
  #negative outlier scores don't make sense for this best score calculation if the direction is positive (generated when outlier frequency is higher in control group; there are a few instances where both the positive outliers and negative outliers are negative >> set to 0 for best score)
    OP_output[["DEscores"]]$best[!is.na( OP_output[["DEscores"]]$best) & (OP_output[["DEscores"]]$NegativeOutlierScores < 0) & (OP_output[["DEscores"]]$PositiveOutlierScores < 0)] <- 0  
  }
  return(OP_output)
}
