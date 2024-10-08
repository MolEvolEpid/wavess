#' Calculate tree summary statistics
#'
#' See `vignette('analyze_output')` for more details.
#'
#' @param tr Phylogeny
#' @param timepoints Vector of time points named by tree tip label. All tip
#' lables must be included in this vector. If you want to exclude certain tips,
#' you must drop them from the tree prior to using this function.
#'
#' @return Tibble including 3 tree summary statistics:
#' - Sackin index ([treebalance::sackinI()])
#' - Mean internal over external branch length ratio ([calc_int_over_ext()])
#' - Timepoint parsimony score ([calc_parsimony()])
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- sample(3, 100, replace = TRUE)
#' names(times) <- tr$tip.label
#' calc_tr_summary_stats(tr, times)
calc_tr_summary_stats <- function(tr, timepoints) {
  check_is_phylo(tr, "tr", rooted = TRUE)
  tibble::tibble(
    stat_name = c("sackin", "int_over_ext", "parsimony_score"),
    stat_value = c(
      treebalance::sackinI(tr),
      calc_int_over_ext(tr),
      calc_parsimony(tr, timepoints)
    )
  )
}

#' Calculate mean internal/external branch length ratio
#'
#' @inheritParams calc_tr_summary_stats
#'
#' @return Mean internal/external branch length ratio
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' calc_int_over_ext(tr)
calc_int_over_ext <- function(tr) {
  check_is_phylo(tr, "tr")
  int_bl <- mean(tr$edge.length[tr$edge[, 2] > ape::Ntip(tr)])
  ext_bl <- mean(tr$edge.length[tr$edge[, 2] <= ape::Ntip(tr)])
  int_bl / ext_bl
}

#' Calculate tree parsimony score for sampling times
#'
#' @inheritParams calc_tr_summary_stats
#'
#' @return Parsimony score of the tree based on the timepoints
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- sample(3, 100, replace = TRUE)
#' names(times) <- tr$tip.label
#' calc_parsimony(tr, times)
calc_parsimony <- function(tr, timepoints) {
  check_is_phylo(tr, "tr", rooted = TRUE)
  if (is.null(names(timepoints)) | !all(names(timepoints) %in% tr$tip.label)) {
    stop(
      "timepoints must be a vector named by tr tip labels, ",
      "and must contain all tip labels."
    )
  }
  phangorn::parsimony(
    tr,
    phangorn::phyDat(factor(timepoints), type = "USER")
  )
}
