#' Calculate tree summary statistics
#'
#' See `vignette('analyze_output')` for more details.
#'
#' @param tr Rooted phylogeny
#' @param timepoints Factor vector of time points named by tree tip label. The
#'   levels should be ordered correctly (most often by sampling time). All tip
#'   lables must be included in this vector. If you want to exclude certain
#'   tips, you must drop them from the tree prior to using this function.
#' @param bl_thresh Branch length threshold under which branches are collapsed
#'   to prior to calculations (default: 1e-08). This is the `tol` argument in
#'   [ape::di2multi()].
#'
#' @return Tibble including 3 tree summary statistics:
#' - Mean leaf depth (normalized Sackin index)
#' - Mean branch length
#' - Mean internal branch length
#' - Mean external branch length
#' - Mean divergence (mean per-generation root-to-tip distance)
#' - Mean diversity (mean per-generation tip-to-tip distance)
#' - Divergence slope across timepoints
#' - Diversity slope across timepoints
#' - Timepoint transition score normalized by the number of timepoints
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- factor(sample(3, 100, replace = TRUE), levels = 1:3)
#' names(times) <- tr$tip.label
#' calc_tr_stats(tr, times)
calc_tr_stats <- function(tr, timepoints, bl_thresh = 1e-08, resolve_timepoints = TRUE) {
  check_is_phylo(tr, "tr")
  if (is.null(names(timepoints)) || !all(names(timepoints) %in% tr$tip.label) || !is.factor(timepoints)) {
    stop(
      "timepoints must be a factor named by tr tip labels, ",
      "and must contain all tip labels."
    )
  }

  tr_poly <- ape::di2multi(tr, tol = bl_thresh)

  transitions <- NA # only 1 timepoint
  if (dplyr::n_distinct(timepoints) > 1) {
    tr_resolved <- tr
    if(!ape::is.binary(tr_poly) & resolve_timepoints){
      # change timepoints to be factors from 1-n
      trait <- factor(as.numeric(timepoints), levels = sort(unique(as.numeric(timepoints))))
      names(trait) <- names(timepoints)
      tr_poly_ordered <- paleotree::resolveTreeChar(tr_poly, trait, orderedChar = TRUE, stateBias = 'primitive')
      tr_resolved <- ape::multi2di(tr_poly_ordered, random = FALSE)
    }
    transitions <- phangorn::parsimony(tr_resolved, phangorn::phyDat(factor(timepoints), type = "USER")) / (dplyr::n_distinct(timepoints) - 1)
  }

  diverg_divers <- lapply(unique(timepoints), function(x) {
    calc_tr_dists(tr, names(timepoints)[timepoints == x]) |>
      dplyr::mutate(timepoint = x, .before = 1)
  }) |>
    dplyr::bind_rows() |>
    dplyr::summarize(
      mean_divergence = mean(.data$mean_rtt, na.rm = TRUE),
      mean_diversity = mean(.data$mean_ttt, na.rm = TRUE)
    ) |>
    unlist() |>
    unname()

  dists <- lapply(unique(timepoints), function(y) {
    tips <- names(timepoints)[timepoints == y]
    ntip <- length(tr$tip.label)
    root_node <- ntip + 1
    d <- ape::dist.nodes(tr)
    colnames(d) <- rownames(d) <- c(tr$tip.label, paste0("node", root_node:(ntip + ape::Nnode(tr))))
    # root-to-tip distance
    rtt <- d[tips, root_node, drop = FALSE] |>
      tibble::enframe() |>
      dplyr::mutate(stat_name = "root_to_tip")
    # tip-to-tip distance
    d <- d[tips, tips]
    ttt <- d[upper.tri(d)] |>
      tibble::enframe() |>
      dplyr::mutate(stat_name = "tip_to_tip") |>
      dplyr::mutate(name = as.character(.data$name))
    if (nrow(ttt) == 0) {
      warning("Generation ", y, " has only one tip, cannot calculate diversity.")
    }
    dplyr::bind_rows(rtt, ttt) |>
      dplyr::mutate(timepoint = y, .before = 1)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(timepoint = as.numeric(as.character(.data$timepoint)))

  slopes <- dplyr::bind_rows(
    stats::lm(value ~ timepoint, dists |>
      dplyr::filter(.data$stat_name == "root_to_tip")) |>
      stats::coef() |>
      tibble::enframe() |>
      dplyr::mutate(type = "root_to_tip"),
    stats::lm(value ~ timepoint, dists |>
      dplyr::filter(.data$stat_name == "tip_to_tip")) |>
      stats::coef() |>
      tibble::enframe() |>
      dplyr::mutate(type = "tip_to_tip")
  ) |>
    dplyr::filter(.data$name == "timepoint")

  tibble::tibble(
    stat_name = c(
      "mean_leaf_depth", "mean_bl",
      "mean_int_bl", "mean_ext_bl",
      "mean_divergence", "mean_diversity",
      "divergence_slope", "diversity_slope",
      "transition_score"
    ),
    stat_value = c(
      treebalance::avgLeafDepI(tr_poly), # average leaf depth (normalized sackin)
      mean(tr$edge.length), # branch lengths
      mean(tr$edge.length[tr$edge[, 2] > ape::Ntip(tr)]), # internal branch lengths
      mean(tr$edge.length[tr$edge[, 2] <= ape::Ntip(tr)]), # external branch lengths
      diverg_divers[1], # divergence
      diverg_divers[2], # diversity
      slopes$value[slopes$type == "root_to_tip"],
      slopes$value[slopes$type == "tip_to_tip"],
      transitions # transition score
    )
  )
}


#' Calculate tree distances (root-to-tip and tip-to-tip)
#'
#' @param tips Tips use when calculating mean distance
#' @inheritParams calc_tr_stats
#'
#' @return Mean root-to-tip distance and mean tip-to-tip distance
#' @noRd
calc_tr_dists <- function(tr, tips = NULL) {
  check_is_phylo(tr, "tr")
  if (is.null(tips)) {
    tips <- tr$tip.label
  }
  ntip <- length(tr$tip.label)
  root_node <- ntip + 1
  d <- ape::dist.nodes(tr)
  colnames(d) <- rownames(d) <- c(tr$tip.label, paste0("node", root_node:(ntip + ape::Nnode(tr))))
  # mean root-to-tip distance
  mean_rtt <- mean(d[tips, root_node])
  # mean tip-to-tip distance
  d <- d[tips, tips]
  mean_ttt <- mean(d[upper.tri(d)])
  return(tibble::tibble(mean_rtt, mean_ttt))
}
