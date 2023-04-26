#' INKA
#'
#' @description
#' Calculates kinase activity score as described in INKA. Combines kinase centric and
#' substrate centric information.
#'
#' @details
#' INKA infers regulator activities
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
#' run_INKA(mat, net, minsize=0)
run_INKA <- function(mat,
                     network,
                     .source = .data$source,
                     .target = .data$target,
                     .mor = .data$mor,
                     minsize = 0
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

    targets <- network_filtered %>%
      dplyr::filter(source == kinase) %>%
      dplyr::pull(target)

    mor_targets <- mat %>%
      dplyr::filter(rownames(mat) %in% targets) %>%
      rownames_to_column("target") %>%
      left_join(network_filtered %>%
                  dplyr::filter(source == kinase), by = "target") %>%
      pull(mor)

    mat_mor <- (mat %>%
                  dplyr::filter(rownames(mat) %in% network_filtered$target[network_filtered$source == kinase])) * mor_targets

  substrate.centric <- mat_mor %>%
      base::colSums()

  mapping_kinase <- read_tsv("resources/kinase_mapping.tsv", col_types = cols())
  kin_id <- mapping_kinase %>%
    filter(external_gene_name == kinase) %>%
    pull(ensembl_gene_id)

  targets.k <- map(kin_id, function(id){
    network_filtered %>%
      dplyr::filter(str_detect(string = target, pattern = id)) %>%
      dplyr::pull(target)

  }) %>% unlist()

  kinase.centric <- mat %>%
    dplyr::filter(rownames(mat) %in% targets.k) %>%
    base::colSums()

    score <-  sign(substrate.centric) * sqrt(abs(substrate.centric) * abs(kinase.centric))

    if (score == 0){
      score <- NA
    }

    data.frame(source = kinase, condition = colnames(mat), score = c(score, substrate.centric, kinase.centric), method = rep(c("INKA", "INKA_substrate_centric", "INKA_kinase_centric"), each = ncol(mat)) )
  })

  rownames(scores) <- NULL
  return(scores)
}
