#' Check define_growth_curve
#'
#' @inheritParams define_growth_curve
#'
#' @return error if input parameters are incorrect
#'
#' @examples
#' check_define_growth_curve_inputs('logistic', 3000, 2000, 1, 30, 0.9)
check_define_growth_curve_inputs <- function(curve_type, gN, K, n0, gS, pK){
  check_is_pos(gN, 'gN')
  check_is_pos(K, 'K')
  check_is_numeric(n0, 'n0')
  check_is_numeric(gS, 'gS')
  check_is_numeric(pK, 'pK')
  if(!curve_type %in% c('logistic', 'linear', 'constant')){
    stop('`curve_type` must be logistic, linear, or constant. You provided: ',
         curve_type)
  }else if(n0 > K){
    stop('n0 must be a number \u2264K, but is ', n0)
  }else if(gS > gN){
    stop('gS must be a number \u2264gN, but is ', gS)
  }else if(pK < 0 | pK > 1){
    stop('pK must be in the range [0,1], but is ', pK)
  }
}

#' Check define_sampling_scheme
#'
#' @inheritParams define_sampling_scheme
#'
#' @return error if input parameters are incorrect
#'
#' @examples
#' check_define_sampling_scheme_inputs(define_growth_curve(), 300, 20)
check_define_sampling_scheme_inputs <- function(growth_curve, sampling_frequency, max_samp){
  check_is_pos(sampling_frequency, 'sampling_frequency')
  check_is_pos(max_samp, 'max_samp')
  check_is_df(growth_curve, 'growth_curve')
  if(!all(colnames(growth_curve) %in% c('generation', 'active_cell_count'))){
    stop('`growth_curve` must contain the columns `generation` and `active_cell_count`, but contains instead: ',
         paste(colnames(growth_curve), collapse = ', '))
  }else if(sampling_frequency > max(growth_curve$generation)){
    stop('sampling_frequency must be a number \u2264maximum generation, but is ', sampling_frequency)
  }
}

#' Check if a variable is numeric
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not numeric
#'
#' @examples
#' check_is_numeric(1, 'test_var')
check_is_numeric <- function(x, var_name){
  if(!is.numeric(x)){
    stop(var_name, ' must be numeric, but is a ', class(x))
  }
}

#' Check if a variable is positive
#'
#' @inheritParams check_is_numeric
#'
#' @return error if variable is not positive
#'
#' @examples
#' check_is_pos(1, 'test_var')
check_is_pos <- function(x, var_name){
  check_is_numeric(x, var_name)
  if(x <= 0){
    stop(var_name, ' must be a positive number, but is ', x)
  }
}

#' Check if a variable is a dataframe or tibble
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not numeric
#'
#' @examples
#' check_is_df(define_growth_curve(), 'test_var')
check_is_df <- function(x, var_name){
  if(!is.data.frame(x)){
    stop(var_name, ' must be a data frame or tibble, but is a ', class(x))
  }
}
