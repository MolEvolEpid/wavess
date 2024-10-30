#' Calculate nucleotide substitution probabilities
#'
#' For this, you should probably use a set of within-host sequences for the
#' genome region of interest.
#'
#' @param aln An alignment of class ape::DNAbin
#' @param tr Starting tree
#' @param model A string providing model (e.g. "GTR+G(4)+I")
#' @param rearrangement Type of tree rearrangements to perform, one of
#' "none", "NNI", "stochastic" or "ratchet" (default: "none")
#'
#' @return
#' Matrix of nucleotide substitution probabilities.
#' Columns are from and rows are to.
#' @export
#'
#' @examples
#' # NOTE: This is just an example.
#' # (**TODO: ADD AN ACTUAL WITHIN-HOST HIV ENV ALIGNMENT AS AN EXAMPLE?**)
#' estimate_q(hiv_env_flt_2022)
estimate_q <- function(aln, tr = NULL, model = "GTR+R(4)+I",
                       rearrangement = "none") {
  check_estimate_q_inputs(aln, tr, model, rearrangement)
  if (is.null(tr)) {
    tr <- ape::bionj(ape::dist.dna(aln, model = "TN93"))
  }
  fit <- phangorn::pml_bb(aln, model, start = tr, rearrangement = rearrangement)
  # it appears that fit$Q isn't only the rates because it's symmetric,
  # so have to multiply by the base frequencies to get the q matrix
  rate_mat <- matrix(c(
    0, fit$Q[1:3], fit$Q[1], 0, fit$Q[4:5], fit$Q[2], fit$Q[4],
    0, fit$Q[6], fit$Q[3], fit$Q[5:6], 0
  ), nrow = 4, ncol = 4)
  bf_mat <- matrix(c(
    fit$bf[1], rep(0, 4), fit$bf[2], rep(0, 4), fit$bf[3],
    rep(0, 4), fit$bf[3]
  ), nrow = 4, ncol = 4)
  q <- rate_mat %*% bf_mat
  rownames(q) <- c("A", "C", "G", "T")
  colnames(q) <- c("A", "C", "G", "T")
  diag(q) <- -rowSums(q)
  # rate_mat <- rate_mat |>
  #   tibble::as_tibble(rownames = "nt_from")
  # }
  return(q)
}

#' Convert q matrix to nucleotide substitution probabilities
#'
#' Probabilities are converted using the equation p = exp(t*Q) For the GTR
#' model, Q should be the transition rate matrix including nucleotide
#' substitution probabilities so that the nucleotide frequencies are maintained
#' over time The diagonal is such that the rows sum to 1. The mutation rate
#' ("t", i.e. `mut_rate`) = 1 means 1 substitution per site. So you want
#' `mut_rate` to be the number of substitutions per site for 1 generation.
#'
#' @param q Q matrix where from nucleotide is rows and to is columns
#' @param mut_rate per-generation overall mutation rate
#'
#' @return
#' Matrix of nucleotide substitution probabilities.
#' Rows are from, columns are to.
#' @noRd
calc_nt_sub_probs_from_q <- function(q, mut_rate) {
  prob_mat <- ape::matexpo(mut_rate * q)
  diag(prob_mat) <- 0
  prob_mat <- prob_mat / rowSums(prob_mat)
  rownames(prob_mat) <- rownames(q)
  colnames(prob_mat) <- colnames(q)
  return(prob_mat)
}
