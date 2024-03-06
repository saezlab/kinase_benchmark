scale_scores <- function(mat, scaling = "sd"){
  if (scaling == "max"){
    map_dfc(colnames(mat), function(col_i){
      mat[col_i]/max(abs(mat[,col_i]), na.rm = T)
    })
  } else if (scaling == "sd") {
    scale(mat, center = FALSE, scale = TRUE)[,]
  } else {
    mat
  }

}
