plotOutlieR <- function(outlieR_Zscores, outlieRplus=F, input_df, genesAsRownames=T, gene_column=NULL, ref_samples, test_samples, mfg, bin_interval=1, Z_interval=1, scale_fit=1, font_size=1, ylimit="default", test_col="default", ref_col="default"){
  if(genesAsRownames==F){
    rownames(input_df) <- input_df[ , gene_column]
  }
  if(test_col=="default"){
    test_col <- rgb(0,139,0,120, maxColorValue = 255)
  }
  if(ref_col=="default"){
    ref_col <- rgb(210,105,30,160, maxColorValue = 255)
  }
  if(genesAsRownames==F){
    rownames(input_df) <- input_df[ , gene_column]
  }
  if(outlieRplus==F){
    outlieR_df <- outlieR_Zscores
  } else {
    outlieR_df <- outlieR_Zscores[[1]]
  }
  outlieR_df <- outlieR_df[mfg, intersect(colnames(outlieR_df),test_samples), drop=F]
  ref_samples <- intersect(ref_samples, colnames(input_df))
  test_samples <- intersect(test_samples, colnames(input_df))
  ref_df <- as.matrix(input_df[mfg, ref_samples, drop=F])
  distr_min <- min(as.numeric(input_df[mfg, c(test_samples, ref_samples)]), na.rm = T)
  distr_max <- max(as.numeric(input_df[mfg, c(test_samples, ref_samples)]), na.rm = T)
  distr_diff <- distr_max - distr_min
  breaks <- seq(floor(distr_min), ceiling(distr_max), by=bin_interval)
  if(ylimit=="default"){
    x = hist(as.matrix(ref_df[mfg, ]), col=ref_col, xlim = c(distr_min - 0.2*distr_diff, distr_max + 0.2*distr_diff), xlab = mfg, breaks = breaks, main = "", cex.axis=font_size, cex.lab=font_size)
  } else {
    x = hist(as.matrix(ref_df[mfg, ]), col=ref_col, xlim = c(distr_min - 0.2*distr_diff, distr_max + 0.2*distr_diff), xlab = mfg, breaks = breaks, main = "", cex.axis=font_size, cex.lab=font_size, ylim = c(0, ylimit))
  }
  xfit <- seq(distr_min - 0.2*distr_diff, distr_max + 0.2*distr_diff, length=50)
  ref_mean <- mean(as.numeric(ref_df[mfg, ]))
  ref_sd <- sd(as.numeric(ref_df[mfg, ]))
  yfit <- dnorm(xfit, mean = ref_mean, sd = ref_sd)
  yfit <- scale_fit * yfit * diff(x$mids[1:2]) * length(ref_samples)
  lines(xfit, yfit, col="black", lwd=2)
  if(ylimit=="default"){
    hist(as.matrix(input_df[mfg, test_samples]), col=test_col, xlim = c(distr_min - 0.2*distr_diff, distr_max + 0.2*distr_diff), breaks=breaks, add=T, cex.axis=font_size, cex.lab=font_size)
  } else {
    hist(as.matrix(input_df[mfg, test_samples]), col=test_col, xlim = c(distr_min - 0.2*distr_diff, distr_max + 0.2*distr_diff), breaks=breaks, add=T, cex.axis=font_size, cex.lab=font_size, ylim = c(0, ylimit))
  }
  Z_min <- round(min(outlieR_df[mfg, ]))
  Z_max <- round(max(outlieR_df[mfg, ]))
  Z_set <- seq(Z_min, Z_max, by=Z_interval)
  for(i in 1:length(Z_set)){
    abline(v=(ref_mean+ref_sd*Z_set[i]), col="red", lwd=2.5)
    text(ref_mean+ref_sd*Z_set[i]-0.2*ref_sd, max(x$counts)-0.3, as.character(Z_set[i]), col = "red", cex = font_size)
  }
  text(distr_min - 0.15*distr_diff, max(x$counts)-0.1, "Z-score", col="red", cex = 1.2 * font_size)
}
