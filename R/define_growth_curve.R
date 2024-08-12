# Get logistic growth curve for infected active cells

#' Define growth curve
#'
#' @param curve_type type of growth, one of:
#' `constant_pop` (constant population size),
#' `linear` (constant growth),
#' `logistic` (logistic growth; default)
#' @param gN final sampling generation (default: 3000)
#' @param K carrying capacity (default: 2000)
#' @param n0 starting population size (default: 1; used for linear and logistic)
#' @param gS seroconversion generation (default: 30; used for logistic)
#' @param pK proportion of carrying capacity population size at seroconversion
#' (default: 0.9; used for logistic)
#'
#' @return tibble with two columns: generation and active cell count
#' @export
#'
#' @examples
#' define_growth_curve()
#' define_growth_curve(curve_type = 'constant')
#' define_growth_curve(gN = 1000)
define_growth_curve <- function(curve_type = 'logistic',
                                gN = 3000,
                                K = 2000,
                                n0 = 1,
                                gS = 30,
                                pK = 0.9){
  check_define_growth_curve_inputs(curve_type, gN, K, n0, gS, pK)
  if(curve_type == 'logistic'){
    gen_df <- tibble::tibble(generation = 0:gN,
                   # infected cell population size over time
                   active_cell_count = get_logistic_curve(n0, K, gS, pK, gN))
  }else if(curve_type == 'constant'){
    gen_df <- tibble::tibble(generation = 0:gN,
                             # infected cell population size over time
                             active_cell_count = K)
  }else if(curve_type == 'linear'){
    gen_df <- tibble::tibble(generation = 0:gN,
                             # infected cell population size over time
                             active_cell_count = floor(seq(n0, K, length.out = gN+1)))
  }
  return(gen_df)
}



#' Get logistic growth curve
#'
#' Note: to get the equations used in this function, I solved for the
#' growth rate and midpoint based on n0, K, gS, and pK.
#'
#' @inheritParams define_growth_curve
#'
#' @return Population size at each generation
get_logistic_curve <- function(n0, K, gS, pK, gN){
  # number of cells at seroconversion generation
  nS <- K*pK
  # useful constants
  c0 <- log(K/n0 - 1)
  cS <- log(K/nS - 1)
  # get growth rate
  k = -(cS-c0)/gS
  # get midpoint
  xM = c0/k
  return(floor(sapply(0:gN, function(x) logistic_fn(x, K, k, xM))))
}


#' Get y value of logistic curve
#'
#' Get the y value for a logistic function given the
#' x value, carrying capacity, growth rate, and midpoint
#'
#' @param x x value of logistic function (generation)
#' @param growth_rate growth rate
#' @param midpoint midpoint of logistic curve
#'
#' @return y value of logistic curve
#'
#' @inheritParams define_growth_curve
logistic_fn <- function(x, K, growth_rate, midpoint){
  K/(1+exp(-growth_rate*(x-midpoint)))
}
