#' Over Representation Analysis (chisq)
#'
#' @description
#' Calculates regulatory activities using chisq.
#'
#' @details
#' chisq measures the overlap between the target feature set and a list of most
#' altered molecular features in mat. The most altered molecular features can
#' be selected from the top and or bottom of the molecular readout distribution,
#' by default it is the top 5% positive values. With these, a contingency table
#' is build and a one-tailed Fisher’s exact test is computed to determine if a
#' regulator’s set of features are over-represented in the selected features
#' from the data. The resulting score, `chisq`, is the minus log10 of the
#' obtained p-value.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param n_up Integer indicating the number of top targets to slice from mat.
#' @param n_bottom Integer indicating the number of bottom targets to slice from
#'  mat.
#' @param n_background Integer indicating the background size of the sliced
#'  targets. If not specified the number of background targets is determined by
#'  the total number of unique targets in the union of `mat` and `network`.
#' @param with_ties Should ties be kept together? The default, `TRUE`,
#'  may return more rows than you request. Use `FALSE` to ignore ties,
#'   and return the first `n` rows.
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
#' @param minsize Integer indicating the minimum number of targets per source.
#' @inheritDotParams stats::fisher.test -x -y
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_chisq(mat, net, minsize=0)
run_chisq <- function(mat,
                    network,
                    .source = source,
                    .target = target,
                    n_up = ceiling(0.05 * nrow(mat)),
                    n_bottom = 0,
                    n_background = 20000,
                    with_ties = TRUE,
                    seed = 42,
                    minsize = 5,
                    ...) {

  # NSE vs. R CMD check workaround
  condition <- p <- p_value <- rand <- score <- source <- statistic <- target    <- targets <- value <- NULL


  # Before to start ---------------------------------------------------------
  network <- network %>%
    rename_net({{ .source }}, {{ .target }})
  network <- filt_minsize(rownames(mat), network, minsize)
  regulons <- extract_sets(network)

  ns <- .chisq_check_ns(n_up, n_bottom, n_background, network, mat)
  n_up <- ns[1]
  n_bottom <- ns[2]
  n_background <- ns[3]

  withr::with_seed(seed, {
    targets <- .chisq_slice_targets(mat, n_up, n_bottom, with_ties)
  })

  # Run analysis ------------------------------------------------------------
  .chisq_analysis(regulons, targets, n_background, ...)
}

# Helper functions --------------------------------------------------------
#' Wrapper to execute `run_chisq()` logic one finished preprocessing of data
#'
#' @inheritParams run_chisq
#' @param regulons Named list; names from `source` and values
#'  from `target`.
#' @param targets Named list; names from columns of `mat` and
#'  values from sliced data of `mat`.
#'
#' @inherit run_scira return
#' @keywords internal
#' @noRd
.chisq_analysis <- function(regulons, targets, n_background, ...) {

  # NSE vs. R CMD check workaround
  p.value <- NULL

  expand_grid(source = names(regulons), condition = names(targets)) %>%
    rowwise(source, condition) %>%
    summarise(.chisq_fisher_exact_test(
      expected = regulons[[source]],
      observed = targets[[condition]],
      n_background = n_background,
      ...
    ),
    .groups = "drop"
    ) %>%
    select(source, condition,
           p_value = p.value
    ) %>%
    mutate(score = -log10(p_value)) %>%
    add_column(statistic = "chisq", .before = 1) %>%
    select(statistic, source, condition, score, p_value)
}

#' Fisher Exact Test
#'
#' @inheritParams run_chisq
#' @inheritParams .chisq_contigency_table
#'
#' @return Single row summary "glance" of a object of class `htest`.
#' @keywords internal
#' @noRd
.chisq_fisher_exact_test <- function(expected, observed, n_background, ...) {
  exec(
    .fn = stats::chisq.test,
    x = .chisq_contingency_table(expected, observed, n_background),
    y = NULL,
    !!!list(...)
  ) %>%
    broom::glance()
}

#' Create contingency table
#'
#' @inheritParams run_chisq
#' @param expected Vector with expected targets
#' @param observed Vector with observed targets
#'
#' @return 2 x 2 matrix
#' @keywords internal
#' @noRd
.chisq_contingency_table <- function(expected, observed, n_background) {
  true_positive <- intersect(observed, expected) %>% length()
  false_positive <- setdiff(expected, observed) %>% length()
  false_negative <- setdiff(observed, expected) %>% length()
  true_negative <- (n_background -
                      true_positive - false_positive - false_negative)

  c(true_positive, false_positive, false_negative, true_negative) %>%
    matrix(nrow = 2, ncol = 2, byrow = FALSE)
}

#' Slice targets per condition
#'
#' @inheritParams run_chisq
#' @return Named list with sliced targets per condition.
#'
#' @keywords internal
#' @noRd
.chisq_slice_targets <- function(mat, n_up, n_bottom, with_ties) {

  # NSE vs. R CMD check workaround
  rand <- targets <- target <- condition <- value <- NULL

  mat %>%
    as_tibble(rownames = "target") %>%
    tidyr::pivot_longer(
      cols = -target,
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(rand=stats::rnorm(n())) %>%
    arrange(condition, value, rand) %>%
    group_by(condition) %>%
    dplyr::do(bind_rows(utils::head(., n = n_bottom), utils::tail(., n = n_up))) %>%
    arrange(condition) %>%
    summarise(
      targets = rlang::set_names(list(target), condition[1]),
      .groups = "drop"
    ) %>%
    pull(targets)
}

#' Check values of variables with n_prefix
#'
#' Set convenient default values for the ns so that downstream
#' functions work fine.
#'
#' @inheritParams run_chisq
#'
#' @return ns modified if necessary.
#'
#' @keywords internal
#' @noRd
.chisq_check_ns <- function(n_up, n_bottom, n_background, network, mat) {
  if (is.null(n_background)) {
    n_background <- network %>%
      pull(target) %>%
      unique() %>%
      union(rownames(mat)) %>%
      length()
  } else if (n_background < 0) {
    abort("`n` must be a non-missing positive number.")
  }

  if (n_up + n_bottom >= nrow(mat)) {
    n_up <- nrow(mat)
    n_bottom <- 0
  }

  c(n_up, n_bottom, n_background)
}
