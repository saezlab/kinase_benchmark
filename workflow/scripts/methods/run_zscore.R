#' RoKAI_zscore
#'
#' @description
#' Calculates Z-score for kinase activity  described in RoKAI.
#'
#' @details
#' RoKAI_zscore infers regulator activities
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `source`: Source nodes of `network`.
#'  2. `condition`: Condition representing each column of `mat`.
#'  3. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_zscore_RoKAI(mat, net, minsize=0)
run_zscore_RoKAI <- function(mat,
                             network,
                             .source = .data$source,
                             .target = .data$target,
                             .mor = .data$mor,
                             minsize = 5
) {
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  network <- network %>%
    rename_net({{ .source }}, {{ .target }}, {{ .mor }})

  network_filtered <- network %>%
    dplyr::filter(target %in% rownames(mat)) %>%
    dplyr::group_by(source) %>%
    dplyr::filter(dplyr::n() >= minsize)

  # Analysis ----------------------------------------------------------------
  kin_sub <- network_filtered %>%
    pivot_wider(values_from = mor, names_from = target, values_fill = 0) %>%
    column_to_rownames("source")

  scores <- purrr::map_dfr(1:ncol(mat), function(i_mat){
    V <- mat[i_mat] %>%
      drop_na()

    S <- sd(V[,1])

    valid_pps <- intersect(rownames(V), colnames(kin_sub))

    kin_sub_f <- kin_sub[, valid_pps]
    V_f <- V %>%
      filter(rownames(V) %in% valid_pps)
    V_f <- V_f[valid_pps,] %>%
      as.matrix()

    kinaseScores <- (as.matrix(kin_sub_f) %*% V_f) / (S * sqrt(abs(as.matrix(kin_sub_f)) %*% rep(1, length(V_f))))
    kinaseScores <- kinaseScores[!is.na(kinaseScores),]

    data.frame(source = names(kinaseScores), condition = colnames(V), score = as.numeric(as.vector(kinaseScores)), method ="RoKAI_z")
  })

  rownames(scores) <- NULL
  return(scores)

}



#' KSEA_zscore
#'
#' @description
#' Calculates Z-score for kinase activity described in KSEA.
#'
#' @details
#' KSEA_zscore infers regulator activities
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `source`: Source nodes of `network`.
#'  2. `condition`: Condition representing each column of `mat`.
#'  3. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_zscore_KSEA(mat, net, minsize=0)
run_zscore_KSEA <- function(mat,
                             network,
                             .source = .data$source,
                             .target = .data$target,
                             .mor = .data$mor,
                             minsize = 5
) {
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  network <- network %>%
    rename_net({{ .source }}, {{ .target }}, {{ .mor }})

  network_filtered <- network %>%
    dplyr::filter(target %in% rownames(mat)) %>%
    dplyr::group_by(source) %>%
    dplyr::filter(dplyr::n() >= minsize)
  network_filtered <- network_filtered[c("source", "target", "mor")] #remove likelihood column

  # Analysis ----------------------------------------------------------------
  # Code is taken from KSEAapp::KSEA.Scores and adjusted to include mode of regulation
  # and different structure of data inputs
  scores <- purrr::map_dfr(1:ncol(mat), function(i_mat){
    mat_c <- mat[i_mat] %>%
      drop_na() %>%
      rownames_to_column("target")

    KSdata <- full_join(network_filtered, mat_c, by = "target") %>%
      drop_na()
    colnames(KSdata)[4] <- "pps" #careful! if data frame structure changes renaming won't work

    KSdata <- KSdata %>%
      dplyr::mutate(value = mor*pps)

    kinase.list <- as.vector(KSdata$source)
    kinase.list <- as.matrix(table(kinase.list))
    Mean.FC <- aggregate(value ~ source, data = KSdata,
                        FUN = mean)
    Mean.FC <- Mean.FC[order(Mean.FC[, 1]), ]
    Mean.FC$mS <- Mean.FC[, 2]
    Mean.FC$Enrichment <- Mean.FC$mS/abs(mean(mat_c[,2], na.rm = T))
    Mean.FC$m <- kinase.list
    Mean.FC$z.score <- ((Mean.FC$mS - mean(mat_c[,2], na.rm = T)) *
                         sqrt(Mean.FC$m))/sd(mat_c[,2], na.rm = T)
    Mean.FC$p.value <- pnorm(-abs(Mean.FC$z.score))
    Mean.FC$FDR <- p.adjust(Mean.FC$p.value, method = "fdr")
    Mean.FC <- Mean.FC[order(Mean.FC$source), -2]

    data.frame(source = Mean.FC$source, condition = colnames(mat_c)[2], score = Mean.FC$z.score, method ="KSEA_z")

  })
  rownames(scores) <- NULL
  return(scores)
}
