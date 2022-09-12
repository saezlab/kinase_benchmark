calcJaccard <- function(listNodes){
  # Calculate jaccard index
  dist <- unlist(lapply(combn(listNodes, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))

  # Transform into data.frame
  jaccardIdx <- as_tibble(cbind(t(combn(length(listNodes),2)), dist))
  jaccardIdx <- jaccardIdx %>%
    add_row(V1 = jaccardIdx$V2, V2 = jaccardIdx$V1, dist = jaccardIdx$dist) %>%
    add_row(V1 = 1:length(listNodes), V2 = 1:length(listNodes), dist = 1)

  # Transform into matrix
  jaccardMat <- matrix(NA, nrow = length(listNodes), ncol = length(listNodes))
  for (i in 1:nrow(jaccardIdx)) {
    jaccardMat[jaccardIdx$V1[i], jaccardIdx$V2[i]] <- jaccardIdx$dist[i]
  }

  colnames(jaccardMat) <- names(listNodes)
  rownames(jaccardMat) <- names(listNodes)

  return(jaccardMat)
}

get_mat_plot <- function(mat, main = 'Jaccard index', palette = 'Reds', reverse = F, cellwidth = 22, cellheight = 22, rownames = T, breaks = NA){
  # Get order of clustered methods
  cor_heat <- pheatmap(mat, cluster_rows = T,
                       cluster_cols = T, silent=T)
  idxRows <- cor_heat$tree_row$order
  idxCols <- cor_heat$tree_col$order
  mat <- mat[idxRows,idxCols]

  # Plot
  paletteLength <- 100
  color <- colorRampPalette((RColorBrewer::brewer.pal(n = 7, name =palette)))(paletteLength)
  if (reverse){
    color <- rev(color)
  }

  if (!is.na(breaks)){
    breaks <- c(seq(min(mat), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(mat)/paletteLength, max(mat), length.out=floor(paletteLength/2)))

  }


  mat_heat <- pheatmap(mat, color = color,
                       display_numbers=F, cluster_rows = F, number_color='black', border_color=NA,
                       cluster_cols = T, na_col=NA, cellwidth = cellwidth, cellheight = cellheight,
                       legend=T, breaks=breaks,
                       show_rownames = rownames, show_colnames = T, main=main,
                       silent=T)
  return(as.ggplot(mat_heat))
}
