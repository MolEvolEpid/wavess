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

#' Calculate Q matrix from nucleotide substitution rates
#'
#' Method:
#' - If needed, convert rates per day to rates per generation
#' - Convert rates per generation to probabilities per generation
#' - Make diagonal such that rows of probability matrix p sum to 1
#' - Convert to q matrix, assuming a mutation rate of mut_rate mutations/site/generation,
#' by solving for q in the equation p = exp(q*mut_rate)
#'
#' @param rates 4x4 matrix of individual nucleotide substitution rates for each
#'   potential substitution (per-site per-generation or per-site per-day) with
#'   the following row and column names: A,C,G,T.
#' @param mut_rate Overall per-site per-generation mutation rate.
#' @param generation_time Generation time (if nucleotide substitution rates are
#'   in days rather than generations; default: NULL, assumes that rate matrix is
#'   per-generation)
#'
#' @return Nucletoide substitution rate matrix Q
#' @export
#'
#' @examples
#' calc_q_from_rates(hiv_mut_rates, 2.4e-5, 1.2)
calc_q_from_rates <- function(rates, mut_rate, generation_time = NULL) {
  check_q_rate(rates, "rates")
  check_is_pos(rates, "rates", ok0 = TRUE)
  check_is_pos(mut_rate, "mut_rate")

  rates_per_gen <- rates
  if (!is.null(generation_time)) {
    check_is_pos(generation_time, "generation_time")
    rates_per_gen <- rates * generation_time
  }
  probs_per_gen <- rate_to_probability(rates_per_gen)
  diag(probs_per_gen) <- 1 - rowSums(probs_per_gen)
  hiv_q_mat <- expm::logm(probs_per_gen) / mut_rate
  rownames(hiv_q_mat) <- colnames(hiv_q_mat) <- colnames(rates)
  return(hiv_q_mat)
}

#' Convert a rate to a probability
#'
#' @param rate Rate to be converted
#' @param time Time period (default: 1, e.g. 1 generation)
#'
#' @return Probability that at least one event occurs in the time period,
#'   assuming a Poisson process.
#' @noRd
rate_to_probability <- function(rate, time = 1) {
  check_is_numeric(rate, "rate")
  check_is_pos(time, "time")
  return(1 - exp(-rate * time))
}

