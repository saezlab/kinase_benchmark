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

  # Analysis ----------------------------------------------------------------
  kinases <- network_filtered$source %>%
    base::unique()
  scores <- purrr::map_dfr(kinases, function(kinase){
    S <- network_filtered %>%
      dplyr::filter(source == kinase) %>%
      base::nrow()

    mor_targets <- mat %>%
      dplyr::filter(rownames(mat) %in% network_filtered$target[network_filtered$source == kinase]) %>%
      rownames_to_column("target") %>%
      left_join(network_filtered %>%
                  dplyr::filter(source == kinase), by = "target") %>%
      pull(mor)

    mat_mor <- (mat %>%
                  dplyr::filter(rownames(mat) %in% network_filtered$target[network_filtered$source == kinase])) * mor_targets


    mean_targets <- mat_mor %>%
      base::colMeans()

    mean_all <- mat %>%
      base::colMeans()

    sdv_all <- mat %>%
      summarise_if(is.numeric, sd)

    z.score <- (sqrt(S)/sdv_all) * (mean_targets - mean_all)

    data.frame(source = kinase, condition = colnames(mat), score = as.numeric(as.vector(z.score[1,])), method ="KSEA_z")
  })

  rownames(scores) <- NULL
  return(scores)
}
