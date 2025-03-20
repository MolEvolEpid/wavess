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
#' The [run_wavess()] function takes probabilities,
#' but sometimes in the literature you find rates instead.
#' This function allows you to convert that rate to a probability.
#'
#' @param rate Rate to be converted
#' @param time Time period (default: 1, e.g. 1 generation)
#'
#' @return
#' Probability that at least one event occurs in the time period
rate_to_probability <- function(rate, time = 1) {
  check_is_numeric(rate, "rate")
  check_is_pos(time, "time")
  return(1 - exp(-rate * time))
}
