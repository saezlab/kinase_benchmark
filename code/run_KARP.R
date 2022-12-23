#' KARP (K-Score)
#'
#' @description
#' Calculates kinase activity score described in KARP (K-score).
#'
#' @details
#' KARP infers regulator activities
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
#' run_KARP(mat, net, minsize=0)
run_KARP <- function(mat,
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
    m <- network_filtered %>%
      dplyr::filter(source == kinase) %>%
      base::nrow()
    t <- network %>%
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

    alpha <- mat_mor %>%
      base::colSums()

    beta <- mat %>%
      base::colSums()

    K.score <- alpha/beta * (m/t)^0.5 * 10^6

    data.frame(source = kinase, condition = colnames(mat), score = K.score, method ="KARP")
  })

 rownames(scores) <- NULL
 return(scores)
}
