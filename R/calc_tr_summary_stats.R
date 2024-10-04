#' Calculate tree summary statistics
#'
#' See `vignette('analyze_output')` for more details.
#'
#' @param tr Phylogeny
#' @param timepoints Vector of time points named by tree tip label
#'
#' @return Tibble including 3 tree summary statistics:
#' - Sackin index normalized by the number of tree tips
#' - Mean internal over external branch length ratio
#' - Mean percent of lineages that survived from one time point to the next
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- sample(3, 100, replace = TRUE)
#' names(times) <- tr$tip.label
#' calc_tr_summary_stats(tr, times)
calc_tr_summary_stats <- function(tr, timepoints) {
  check_is_phylo(tr)
  check_is_pos(timepoints, ok0 = TRUE)
  tibble::tibble(
    stat_name = c("sackin_norm", "int_over_ext", "mean_prop_survived"),
    stat_value = c(
      calc_sackin_norm(tr), calc_int_over_ext(tr),
      mean(calc_prop_survived(tr, timepoints)$prop_survived,
        na.rm = TRUE
      )
    )
  )
}

#' Calculate Sackin indexed normalized by number of tips
#'
#' @inheritParams calc_tr_summary_stats
#'
#' @return Normalized Sackin index
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' calc_sackin_norm(tr)
calc_sackin_norm <- function(tr) {
  check_is_phylo(tr)
  treebalance::sackinI(tr) / ape::Ntip(tr)
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
  check_is_phylo(tr)
  int_bl <- mean(tr$edge.length[tr$edge[, 2] > ape::Ntip(tr)])
  ext_bl <- mean(tr$edge.length[tr$edge[, 2] <= ape::Ntip(tr)])
  int_bl / ext_bl
}

#' Calculate percent of lineages that survived from one generation to the next
#'
#' @inheritParams calc_tr_summary_stats
#'
#' @return Tibble including:
#' - `n_seqs`: Number of sequences included in that time point
#' - `n_clusters`: Number of monophyletic clusters for that time point
#' - `prop_survived`: Proportion of lineages that survived from the
#'    previous generation
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' times <- sample(3, 100, replace = TRUE)
#' names(times) <- tr$tip.label
#' calc_prop_survived(tr, times)
calc_prop_survived <- function(tr, timepoints) {
  check_is_phylo(tr)
  check_is_pos(timepoints, ok0 = TRUE)
  storage.mode(timepoints) <- "integer"
  timepts_unique <- unique(timepoints)
  # can't find clusters for first timepoint
  lapply(timepts_unique, function(t) {
    tr_sub <- ape::keep.tip(tr, names(timepoints)[timepoints <= t])
    timepoints_sub <- timepoints[timepoints <= t]
    get_clusters(tr_sub, timepoints_sub, bootstrap = NULL)$pure_subtree_info |>
      dplyr::filter(.data$timepoint == t) |>
      dplyr::group_by(.data$timepoint) |>
      dplyr::summarize(
        n_seqs = sum(.data$subtr_size),
        n_clusters = dplyr::n()
      )
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(prop_survived = .data$n_clusters / dplyr::lag(.data$n_seqs))
}

#' Get monophyletic clusters on the phylogeny
#'
#' This code was modified from https://github.com/Snitkin-Lab-Umich/regentrans/,
#' which is under an MIT license
#'
#' @inheritParams calc_tr_summary_stats
#' @param pureness How pure each cluster should be (must be > 0.5) (optional,
#'   default = 1)
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges
#'   (optional, default = NULL)
#' @param grps Groups from which to keep only one isolate (default: NULL, all
#'   kept)
#'
#' @return list where pure_subtree_info is a data.frame of clusters on
#'   phylogeny, index indicates which element that cluster is in the list of
#'   subtrees, NA indicates it is not part of a subtree; subtrees is an object
#'   of the actual subtrees (can be used for plotting); cluster_pureness is the
#'   purness of each cluster
#' @noRd
get_clusters <- function(tr, timepoints, pureness = 1, bootstrap = NULL,
                         grps = NULL) {
  # get the names of the things in common
  isolates <- intersect(tr$tip.label, names(timepoints))
  if (!is.null(grps)) {
    # subset to one sequence per group
    grp_df <- as.data.frame(grps) |>
      dplyr::distinct(grps, .keep_all = TRUE) |>
      rownames() |>
      as.vector()
    isolates <- intersect(isolates, grp_df)
  }
  # subset timepoints
  timepoints_sub <- timepoints[isolates]
  # subset the tree
  tr <- ape::keep.tip(tr, isolates)

  subtrs_sub <- ape::subtrees(tr)
  pure_subtrees <- get_largest_subtree(
    subtrs = subtrs_sub,
    isolate_labels = timepoints_sub,
    bootstrap = bootstrap,
    pureness = pureness
  )
  pure_subtr_info <- dplyr::bind_cols(
    timepoint = timepoints_sub,
    subtr_size = unlist(pure_subtrees$largest_st),
    index = unlist(pure_subtrees$largest_st_i),
    isolate_id = names(timepoints_sub)
  )
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info |>
    dplyr::mutate(subtr_size = ifelse(.data$subtr_size == 0 & .data$index == 1,
      1, .data$subtr_size
    ))
  # change index from 1 to NA
  pure_subtr_info <- pure_subtr_info |>
    dplyr::mutate(index = ifelse(.data$index == 1, NA, .data$index))
  # add a column to indicate the isolate name if the index = NA
  pure_subtr_info <- pure_subtr_info |>
    dplyr::mutate(isolate_id = ifelse(is.na(.data$index), .data$isolate_id, NA))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info$index) |
    pure_subtr_info$subtr_size == 1, ]


  # potential returns if we want to return subtrees too
  returns <- list(
    pure_subtree_info = pure_subtr_info,
    subtrees = subtrs_sub,
    cluster_pureness = pureness
  )

  return(returns)
}


#' Get largest pure subtrees
#'
#' This code was modified from https://github.com/Snitkin-Lab-Umich/regentrans/,
#' which is under an MIT license
#'
#' @param subtrs Subtrees created using ape::subtrees to look for clustering on.
#'   Should include all isolates of interest.
#' @param isolate_labels Named vector of labels by which pure clusters are
#'   defined. Names must be equivalent to tree tip label names.
#' @param control_labels Named vector of labels known to cluster. Names must be
#'   equivalent to tree tip label names. This controls for clustering by
#'   requiring that the pure clusters must contain multiple of the control
#'   labels.
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges
#'   (keeps > bootstrap; NULL = keep all; default: 90).
#' @param pureness How pure the subtree has to be to call it a "pure" subtree
#'   (default: 1; range 0-1).
#'
#' @return list containing the largest pure subtree that each isolate belongs
#'   to, the index of that subtree, and the edges in that subtree.
#' @noRd
get_largest_subtree <- function(subtrs, isolate_labels, control_labels = NULL,
                                bootstrap = 90, pureness = 1) {
  largest_st_info <- lapply(names(isolate_labels), function(i) {
    # DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE
    # DEFINED AS: 1) HAVE ONLY A SINGLE EPI LABEL, 2) HAVE BOOTSTRAP SUPPORT
    # GREATER THAN bootstrap, 3) INCLUDE MORE THAN ONE CONTROL LABEL
    sts <- sapply(subtrs, FUN = function(st) {
      i_in_subtree <- i %in% st$tip.label # isolate is in subtree
      st_labs <- isolate_labels[intersect(st$tip.label, names(isolate_labels))]
      one_label <- length(unique(st_labs)) == 1 # only one label in subtree
      labs_tab <- table(st_labs)
      labs_max <- max(labs_tab)
      labs_max_lab <- names(labs_tab)[which.max(labs_tab)]
      if (pureness != 1) { # allow some contamination based on pureness
        one_label <- labs_max / length(st_labs) > pureness
      }
      good_bootstrap <- rep(TRUE, length(st$node.label[[1]]))
      if (!is.null(bootstrap)) {
        good_bootstrap <- !is.na(as.numeric(st$node.label[[1]])) &&
          as.numeric(st$node.label[[1]]) > bootstrap
      }
      # alwaystrue if not controlling for another variable
      multiple_control <- ifelse(is.null(control_labels), TRUE,
        length(unique(control_labels[intersect(
          st$tip.label,
          names(control_labels)
        )])) > 1
      ) # more than one control label in subtree
      if (i_in_subtree && one_label && good_bootstrap && multiple_control &&
        isolate_labels[i] == labs_max_lab) {
        length(intersect(
          names(isolate_labels[isolate_labels == labs_max_lab]),
          st$tip.label
        ))
      } else {
        0
      }
    })
    # GET THE LARGEST SUBTREE
    largest_st <- max(sts)
    # GET THE INDEX OF THE LARGEST SUBTREE
    largest_st_i <- which.max(sts)
    # GET EDGES BELONGING TO SUBTREES
    largest_st_edges <- ape::which.edge(
      subtrs[[1]],
      subtrs[[largest_st_i]]$tip.label
    )
    return(list(
      largest_st = largest_st,
      largest_st_i = largest_st_i,
      largest_st_edges = largest_st_edges
    ))
  }) # end sapply
  # reverse lists
  largest_st_info <- reverse_list_str(largest_st_info)
  return(largest_st_info)
} # end get_largest_subtree

#' Reverse list structure
#'
#' This code was modified from https://github.com/Snitkin-Lab-Umich/regentrans/,
#' which is under an MIT license
#'
#' @param ls list you want to reverse
#'
#' @return reversed list
#' @details Reference with example:
#'   https://stackoverflow.com/questions/15263146/revert-list-structure
#' @noRd
reverse_list_str <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list)
}
