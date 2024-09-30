#' Convert a rate to a probability
#'
#' The `wavess()` function takes probabilities,
#' but sometimes in the literature you find rates instead.
#' This function allows you to convert that rate to a probability.
#'
#' @param rate Rate to be converted
#' @param time Time period (default: 1, e.g. 1 generation)
#'
#' @return
#' Probability that at least one event occurs in the time period
#' @export
#'
#' @examples
#' rate_to_probability(2.1)
#' rate_to_probability(0.1)
#' rate_to_probability(3.5e-5)
rate_to_probability <- function(rate, time = 1){
  check_is_numeric(rate, 'rate')
  check_is_pos(time, 'rate')
  return(1-exp(-rate*time))
}
