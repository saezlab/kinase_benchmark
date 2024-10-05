#' RoKAI_linear_model
#'
#' @description
#' Calculates kinase activity with linear model described in RoKAI.
#'
#' @RoKAI_linear_model
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
#' RoKAI_linear_model(mat, net, minsize=0)
run_lm_rokai <- function(mat,
                         network,
                         .source = .data$source,
                         .target = .data$target,
                         .mor = .data$mor,
                         minsize = 5,
                         k = 0.1
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
  # Change format
  mat <- mat %>% filter(rownames(mat) %in% network_filtered$target)

  network_wide <- network_filtered %>%
    dplyr::select(source, target, mor) %>%
    pivot_wider(names_from = target, values_from = mor) %>%
    column_to_rownames("source")

  network_wide[is.na(network_wide)] <- 0
  network_wide <- network_wide[,match(rownames(mat), colnames(network_wide))]

  scores <- map_dfr(1:ncol(mat), function(exp){
    V <- mat[,1] %>% as.matrix()
    rownames(V) <- rownames(mat)

    linear_kinase_activity_inference <- function(V, substrates, k=0.1) {
      validKinases <- rowSums(substrates) > 0
      X <- t(substrates[validKinases, ])
      Y <- V
      n <- nrow(X)
      m <- ncol(X)
      pseudo <- sqrt(k) * diag(m + 1)
      Xplus <- rbind(cbind(rep(1, n), X), pseudo)
      Yplus <- rbind(Y, matrix(rep(0, m + 1)))
      beta <- qr.solve(Xplus, Yplus)
      kinaseBeta <- rep(NA, nrow(substrates))
      kinaseBeta[validKinases] <- beta[-1]
      return(kinaseBeta)
    }

    act <- linear_kinase_activity_inference(V = V, substrates = network_wide, k = k)

    data.frame(source = rownames(network_wide), condition = colnames(mat)[exp], score = act, method ="lmRoKAI")
  })

  rownames(scores) <- NULL
  return(scores)
}
