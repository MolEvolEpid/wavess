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
#' calc_nt_sub_probs(hiv_env_flt_2021)
calc_nt_sub_probs <- function(aln, tr = NULL, model = "GTR+R(4)+I",
                              rearrangement = "none") {
  check_calc_nt_sub_probs_inputs(aln, tr, model, rearrangement)
  if (is.null(tr)) {
    tr <- ape::bionj(ape::dist.dna(aln, model = "TN93"))
  }
  fit <- phangorn::pml_bb(aln, model, start = tr, rearrangement = rearrangement)
  # it appears that fit$Q is actually the rate matrix because it's symmetric,
  # so have to multiply by the base frequencies to get the q matrix
  rate_mat <- matrix(c(
    0, fit$Q[1:3], fit$Q[1], 0, fit$Q[4:5], fit$Q[2], fit$Q[4],
    0, fit$Q[6], fit$Q[3], fit$Q[5:6], 0
  ), nrow = 4, ncol = 4)
  bf_mat <- matrix(c(
    fit$bf[1], rep(0, 4), fit$bf[2], rep(0, 4), fit$bf[3],
    rep(0, 4), fit$bf[3]
  ), nrow = 4, ncol = 4)
  q_mat <- rate_mat %*% bf_mat
  prob_mat <- q_mat / rowSums(q_mat)
  rownames(prob_mat) <- c("A", "C", "G", "T")
  colnames(prob_mat) <- c("A", "C", "G", "T")
  return(data.frame(prob_mat))
}
