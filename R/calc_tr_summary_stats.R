#' Calculate tree summary statistics
#'
#' See `vignette('analyze_output')` for more details.
#'
#' @param tr Rooted phylogeny
#' @param timepoints Vector of time points named by tree tip label. All tip
#' lables must be included in this vector. If you want to exclude certain tips,
#' you must drop them from the tree prior to using this function.
#'
#' @return Tibble including 3 tree summary statistics:
#' - Sackin index ([treebalance::sackinI()])
#' - Mean internal branch length ([calc_int_over_ext()])
#' - Mean external branch length ([calc_int_over_ext()])
#' - Mean internal over external branch length ratio ([calc_int_over_ext()])
#' - Timepoint parsimony score ([calc_parsimony()])
#' - Mean root-to-tip distance ([calc_tr_dists()])
#' - Mean tip-to-tip distance ([calc_tr_dists()])
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- sample(3, 100, replace = TRUE)
#' names(times) <- tr$tip.label
#' calc_tr_summary_stats(tr, times)
calc_tr_summary_stats <- function(tr, timepoints) {
  check_is_phylo(tr, "tr")
  tibble::tibble(
    stat_name = c("sackin", "int_bl", "ext_bl", "int_over_ext", "parsimony_score", "root_to_tip", "tip_to_tip"),
    stat_value = c(
      treebalance::sackinI(tr),
      calc_int_over_ext(tr),
      calc_parsimony(tr, timepoints),
      calc_tr_dists(tr)
    )
  )
}

#' Calculate mean internal/external branch length ratio
#'
#' @inheritParams calc_tr_summary_stats
#'
#' @return Mean internal branch length, mean external branch length, and mean
#'   internal/external branch length ratio
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' calc_int_over_ext(tr)
calc_int_over_ext <- function(tr) {
  check_is_phylo(tr, "tr")
  int_bl <- mean(tr$edge.length[tr$edge[, 2] > ape::Ntip(tr)])
  ext_bl <- mean(tr$edge.length[tr$edge[, 2] <= ape::Ntip(tr)])
  c(int_bl, ext_bl, int_bl / ext_bl)
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
  check_is_phylo(tr, "tr")
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

#' Calculate tree distances (root-to-tip and tip-to-tip)
#'
#' @param tips Tips use when calculating mean distance
#' @inheritParams calc_tr_summary_stats
#'
#' @return Mean root-to-tip distance and mean tip-to-tip distance
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' calc_tr_dists(tr)
calc_tr_dists <- function(tr, tips = NULL){
  check_is_phylo(tr, "tr")
  if(is.null(tips)){
    tips <- tr$tip.label
  }
  ntip <- length(tr$tip.label)
  root_node <- ntip+1
  d <- ape::dist.nodes(tr)
  colnames(d) <- rownames(d) <- c(tr$tip.label, paste0('node', root_node:(ntip+ape::Nnode(tr))))
  # mean root-to-tip distance
  # dv <- mean(d[1:ntip,root_node])
  dv <- mean(d[tips,root_node])
  # mean tip-to-tip distance
  # d <- d[1:ntip,1:ntip]
  d <- d[tips,tips]
  ds <- mean(d[upper.tri(d)])
  return(c(dv, ds))
}
