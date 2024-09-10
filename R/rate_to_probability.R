#' Convert a rate to a Poisson probability
#'
#' The `wavess()` function takes probabilities,
#' but sometimes in the literature you find rates instead.
#' This funciton allows you to convert that rate to a probability.
#'
#' @param rate Rate to be converted
#' @param k Number of events observed in the rate-defined time period (default: 1)
#'
#' @return
#' Probability that the event occurs in the rate-defined time period
#' @export
#'
#' @examples
#' rate_to_probability(2.1)
#' rate_to_probability(0.1)
#' rate_to_probability(3.5e-5)
rate_to_probability <- function(rate, k = 1){
  check_is_numeric(rate, 'rate')
  check_is_numeric(k, 'rate')
  if(k < 0) check_is_pos(k, 'k')
  return(rate^k*exp(-rate)/factorial(k))
}
