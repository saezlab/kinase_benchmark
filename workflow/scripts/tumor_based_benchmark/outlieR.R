outlieR <- function(input_df, ref_samples, test_samples, scaling_sample_set1=NA, scaling_sample_set2=NA, mfg=NA, genesAsRownames=T, gene_column=NULL, z_thresh=2.3265, z_method="normal", scaling_factor=1){
  if(genesAsRownames==F){
    rownames(input_df) <- input_df[ , gene_column]
  }
  ref_samples <- intersect(ref_samples, colnames(input_df))
  if(!is.na(mfg)){
    if(sum(!is.na(input_df[mfg, ref_samples])) > 2){
      input_df <- input_df[mfg, , drop=F]
    } else {
      warning('not enough reference data for mfg')
      return(NA)
    }
  }
  
  #make dataframe containing just values for reference samples
  ref_df <- input_df[ , ref_samples, drop=F]
  
  test_samples <- intersect(test_samples, colnames(input_df))
  test_df <- as.matrix(input_df[ , test_samples, drop=F])
  
  if(is.na(scaling_sample_set1)){
    scaling_samples <- ref_samples
    scaling_df <- ref_df
    scaling_sets <- 1
  } else {
    if(is.na(scaling_sample_set2)){    
      scaling_samples <- intersect(scaling_sample_set1, colnames(input_df))
      scaling_df <- input_df[ , scaling_samples]
      scaling_sets <- 1
    } else {
      scaling_sample_set1 <- intersect(scaling_sample_set1, colnames(input_df))
      scaling_sample_set2 <- intersect(scaling_sample_set2, colnames(input_df))
      scaling_samples <- c(scaling_sample_set1, scaling_sample_set2)
      scaling_df <- input_df[ , scaling_samples]
      scaling_sets <- 2
    }
  }
  
  #need data for at least 3 samples in this dataframe to properly define a distribution; drop genes that don't have at least 3 reference measurements
  ref_df <- ref_df[rowSums(!is.na(ref_df))>2, , drop=F]
  scaling_df <- scaling_df[rownames(ref_df), , drop=F]
  input_df <- input_df[rownames(ref_df), , drop=F]
  test_df <- test_df[rownames(ref_df), , drop=F]
  
  #if method chosen is normal, calculate Z-scores based off of mean and SD and run Wilk-Shapiro test on reference data frame to test if the distribution is actually normal (normal_distribution = 1)
  if(z_method=="normal"){
    ref_df$mean <- rowMeans(ref_df, na.rm = T)
    ref_df <- transform(as.data.frame(ref_df), SD=apply(scaling_df, 1, sd, na.rm = TRUE))
    normality_test <- apply(ref_df[ , 1:(ncol(ref_df)-1)], 1, shapiro.test)
    ref_df$normality_test_pval <- 0
    for(i in 1:nrow(ref_df)){
      ref_df$normality_test_pval[i] <- normality_test[[rownames(ref_df)[i]]]$p.value
    }
    #test_df <- as.matrix(input_df[ , test_samples])
    #rownames(test_df) <- rownames(input_df)
    outlier_Zscores <- (test_df - ref_df$mean)/(ref_df$SD)
    outlier_Zscores <- as.data.frame(outlier_Zscores)
    outlier_Zscores$normal_distribution <- 0
    outlier_Zscores$normal_distribution <- as.numeric(ref_df$normality_test_pval >= 0.05)
    #if method chosen is medMAD, calculate Z-scores based off of median and MAD 
  } else if(z_method=="medMAD"){
    ref_df <- transform(as.data.frame(ref_df), median=apply(ref_df, 1, median, na.rm = TRUE))
    if(scaling_sets==1){
      ref_df <- transform(as.data.frame(ref_df), MAD=apply(scaling_df, 1, mad, na.rm = TRUE))
    } else {
      med_s1 <- apply(scaling_df[ , scaling_sample_set1], 1, median, na.rm = TRUE)
      med_s2 <- apply(scaling_df[ , scaling_sample_set2], 1, median, na.rm = TRUE)
      med_devs1 <- scaling_df[ , scaling_sample_set1] - med_s1
      med_devs2 <- scaling_df[ , scaling_sample_set2] - med_s2
      ref_df <- transform(as.data.frame(ref_df), MAD=apply(cbind(med_devs1, med_devs2), 1, median, na.rm = TRUE))
    }
    #test_df <- as.matrix(input_df[ , intersect(test_samples, colnames(input_df))])
    #rownames(test_df) <- rownames(input_df)
    outlier_Zscores <- (test_df - ref_df$median)/(ref_df$MAD * scaling_factor)
    outlier_Zscores <- as.data.frame(outlier_Zscores)
  }
  
  #determine number of high and low outliers (with Zscores above or below z_thresh, respectively)
  #outlier_Zscores <- as.data.frame(outlier_Zscores)
  outlier_Zscores$Number_of_high_outliers <- rowSums(outlier_Zscores[ , test_samples] > z_thresh, na.rm = T)
  outlier_Zscores$Number_of_low_outliers <- rowSums(outlier_Zscores[ , test_samples] < -z_thresh, na.rm = T)
  
  return(outlier_Zscores)
}
