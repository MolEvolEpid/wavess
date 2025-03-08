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
#' - Mean leaf depth (normalized Sackin index)
#' - Mean internal branch length
#' - Mean external branch length
#' - Timepoint transition score
#' - Mean tip-to-tip distance
#' - Mean divergence (per-generation root-to-tip distance)
#' - Mean diversity (per-generation tip-to-tip distance)
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- sample(3, 100, replace = TRUE)
#' names(times) <- tr$tip.label
#' calc_tr_stats(tr, times)
calc_tr_stats <- function(tr, timepoints) {
  check_is_phylo(tr, "tr")

  if (is.null(names(timepoints)) | !all(names(timepoints) %in% tr$tip.label)) {
    stop(
      "timepoints must be a vector named by tr tip labels, ",
      "and must contain all tip labels."
    )
  }

  transitions <- phangorn::parsimony(tr, phangorn::phyDat(factor(timepoints), type = "USER"))/dplyr::n_distinct(timepoints)

  rtt_ttt <- calc_tr_dists(tr)
  diverg_divers <- lapply(unique(timepoints), function(x) {
    calc_tr_dists(tr, names(timepoints)[timepoints == x]) |>
      dplyr::mutate(timepoint = x, .before = 1)
  }) |>
    dplyr::bind_rows()  |>
    dplyr::summarize(
      mean_divergence = mean(mean_rtt, na.rm = TRUE),
      mean_diversity = mean(mean_ttt, na.rm = TRUE),
      var_divergence = var(mean_rtt, na.rm = TRUE),
      var_diversity = var(mean_ttt, na.rm = TRUE)
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
    rtt <- d[tips, root_node] |>
      tibble::enframe() |>
      dplyr::mutate(stat_name = 'root_to_tip')
    # tip-to-tip distance
    d <- d[tips, tips]
    ttt <- d[upper.tri(d)] |>
      tibble::enframe() |>
      dplyr::mutate(stat_name = 'tip_to_tip') |>
      dplyr::mutate(name = as.character(name))
    dplyr::bind_rows(rtt, ttt) |>
      dplyr::mutate(timepoint = y, .before = 1)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(timepoint = as.numeric(as.character(timepoint)))

  slopes <- bind_rows(lm(value ~ timepoint, dists |>
                 filter(stat_name == 'root_to_tip')) |>
              coef() |>
              tibble::enframe() |>
              dplyr::mutate(type = 'root_to_tip'),
            lm(value ~ timepoint, dists |>
                 dplyr::filter(stat_name == 'tip_to_tip')) |>
              coef() |>
              tibble::enframe() |>
              dplyr::mutate(type = 'tip_to_tip')) |>
    dplyr::filter(name == 'timepoint')

  tibble::tibble(
    stat_name = c("mean_leaf_depth", "corrected_colless", "mean_bl",
                  "mean_int_bl", "mean_ext_bl",
                  "mean_tip_to_tip", "mean_root_to_tip",
                  "mean_divergence", "mean_diversity",
                  "divergence_slope", "diversity_slope",
                  "transition_score"),
    stat_value = c(
      treebalance::avgLeafDepI(tr), # average leaf depth (normalized sackin)
      treebalance::collessI(ape::multi2di(tr), method = 'corrected'), # corrected colless to enable comparison between trees
      mean(tr$edge.length), # branch lengths
      mean(tr$edge.length[tr$edge[, 2] > ape::Ntip(tr)]), # internal branch lengths
      mean(tr$edge.length[tr$edge[, 2] <= ape::Ntip(tr)]), # external branch lengths
      rtt_ttt$mean_ttt, # tip-to-tip
      rtt_ttt$mean_rtt, # root-to-tip
      diverg_divers[1], # divergence
      diverg_divers[2], # diversity
      slopes$value[slopes$type == 'root_to_tip'],
      slopes$value[slopes$type == 'tip_to_tip'],
      transitions # transition score
    )
  )
}
#'
#' #' Calculate mean internal/external branch length ratio
#' #'
#' #' @inheritParams calc_tr_stats
#' #'
#' #' @return Mean internal branch length, mean external branch length, and mean
#' #'   internal/external branch length ratio
#' #' @export
#' #'
#' #' @examples
#' #' tr <- ape::rtree(100)
#' #' calc_int_over_ext(tr)
#' calc_int_over_ext <- function(tr) {
#'   check_is_phylo(tr, "tr")
#'   int_bl <- mean(tr$edge.length[tr$edge[, 2] > ape::Ntip(tr)])
#'   ext_bl <- mean(tr$edge.length[tr$edge[, 2] <= ape::Ntip(tr)])
#'   c(int_bl, ext_bl, int_bl / ext_bl)
#' }
#'
#' #' Calculate tree transition score for sampling times
#' #'
#' #' @inheritParams calc_tr_stats
#' #'
#' #' @return Number of transitions between timepoints on the phylogeny based on
#' #'   parsimony
#' #' @export
#' #'
#' #' @examples
#' #' tr <- ape::rtree(100)
#' #' times <- sample(3, 100, replace = TRUE)
#' #' names(times) <- tr$tip.label
#' #' calc_transition(tr, times)
#' calc_transition <- function(tr, timepoints) {
#'   check_is_phylo(tr, "tr")
#'   if (is.null(names(timepoints)) | !all(names(timepoints) %in% tr$tip.label)) {
#'     stop(
#'       "timepoints must be a vector named by tr tip labels, ",
#'       "and must contain all tip labels."
#'     )
#'   }
#'   phangorn::parsimony(
#'     tr,
#'     phangorn::phyDat(factor(timepoints), type = "USER")
#'   )
#' }

#' Calculate tree distances (root-to-tip and tip-to-tip)
#'
#' @param tips Tips use when calculating mean distance
#' @inheritParams calc_tr_stats
#'
#' @return Mean root-to-tip distance and mean tip-to-tip distance
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
