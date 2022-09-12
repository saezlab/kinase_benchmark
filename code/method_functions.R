runKARP <- function(t.values, PKN, min.targets = 5){
  colnames(t.values) <- c("TARGET", "t")
  t.values <- dplyr::left_join(t.values, PKN, by = "TARGET")
  kinases <- names(table(t.values$GENE))[table(t.values$GENE) >= min.targets] #select kinases with at least n measured targets

  K.scores_df <- map_dfr(kinases, function(kinase){
    m <- t.values %>% dplyr::filter(GENE == kinase) %>% nrow()
    t <- PKN %>% dplyr::filter(GENE == kinase) %>% nrow()

    alpha <- t.values %>%
      dplyr::filter(GENE == kinase) %>%
      dplyr::select(t) %>% sum()

    beta <- t.values %>%
      dplyr::select(t) %>%
      sum()

    K.score <- alpha/beta * (m/t)^0.5 * 10^6

    data.frame(kinase = kinase, K.score = K.score)
  })

  return(K.scores_df %>% dplyr::arrange(dplyr::desc(K.score)))
}
