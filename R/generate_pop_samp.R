#' Generate population growth curve and sampling scheme
#'
#' @param curve_type type of growth, one of:
#' `constant_pop` (constant population size),
#' `logistic` (logistic growth; default)
#' @param gN final sampling generation (default: 3000)
#' @param K carrying capacity (default: 2000)
#' @param n0 starting population size (default: 1; used for logistic)
#' @param g50 generation at which the population size reaches 50% of K
#' (default: 25; used for logistic)
#' @param sampling_frequency frequency in generations at which to record
#' sequences (and counts) (default: 300 generations)
#' @param max_samp maximum number of cells (and thus sequences) to sample in a
#' given generation (default: 20 sequences)
#'
#' @return tibble with two columns:
#' - `generation`: Each generation to be simulated
#' - `active_cell_count`: Number of active cells in each generation
#' - `n_sample_active` Number of sequences to sample in each generation
#' @export
#'
#' @examples
#' generate_pop_samp()
generate_pop_samp <- function(curve_type = "logistic",
                              gN = 3000,
                              K = 2000,
                              n0 = 1,
                              g50 = 25,
                              sampling_frequency = 300,
                              max_samp = 20) {
  check_generate_pop_samp_inputs(curve_type, gN, K, n0, g50, sampling_frequency, max_samp)
  define_growth_curve(curve_type, gN, K, n0, g50, 0.5) |>
    define_sampling_scheme(sampling_frequency, max_samp)
}

# Get logistic growth curve for infected active cells

#' Define growth curve
#'
#' @param gpK generation at which the population size reaches pK of K
#' (default: 25; used for logistic)
#' @param pK proportion of carrying capacity population size at gpK
#' (default: 0.5; used for logistic)
#' @inheritParams generate_pop_samp
#'
#' @return tibble with two columns: generation and active cell count
#' @noRd
define_growth_curve <- function(curve_type = "logistic",
                                gN = 3000,
                                K = 2000,
                                n0 = 1,
                                gpK = 25,
                                pK = 0.5) {
  if (curve_type == "logistic") {
    gen_df <- tibble::tibble(
      generation = 0:gN,
      # infected cell population size over time
      active_cell_count = get_logistic_curve(n0, K, gpK, pK, gN)
    )
  } else if (curve_type == "constant") {
    gen_df <- tibble::tibble(
      generation = 0:gN,
      # infected cell population size over time
      active_cell_count = K
    )
  }
  return(gen_df)
}

#' Get logistic growth curve
#'
#' Note: to get the equations used in this function, I solved for the
#' growth rate and midpoint based on n0, K, gpK, and pK.
#'
#' @param gpK generation at which the population size reaches pK of K
#' (default: 25; used for logistic)
#' @param pK proportion of carrying capacity population size at gpK
#' (default: 0.5; used for logistic)
#' @inheritParams generate_pop_samp
#'
#' @return Population size at each generation
#' @noRd
get_logistic_curve <- function(n0, K, gpK, pK, gN) {
  # number of cells at seroconversion generation
  nS <- K * pK
  # useful constants
  c0 <- log(K / n0 - 1)
  cS <- log(K / nS - 1)
  # get growth rate
  k <- -(cS - c0) / gpK
  # get midpoint
  xM <- c0 / k
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
#' @inheritParams generate_pop_samp
#' @noRd
logistic_fn <- function(x, K, growth_rate, midpoint) {
  K / (1 + exp(-growth_rate * (x - midpoint)))
}

#' Define sampling scheme
#'
#' @param growth_curve output from [define_growth_curve()], or customized
#' growth curve with columns named `generation` and `active_cell_count`
#' @inheritParams generate_pop_samp
#'
#' @return input growth curve tibble with one additional column (`n_sample_active`)
#' containing the number of sequences from active cells to samples
#' @noRd
define_sampling_scheme <- function(growth_curve,
                                   sampling_frequency = 300,
                                   max_samp = 20) {
  ss <- growth_curve |>
    dplyr::rowwise() |>
    # sample up to max_samp sequences every sampling_frequency generations
    dplyr::mutate(n_sample_active = ifelse(.data$generation %% sampling_frequency == 0,
      min(max_samp, .data$active_cell_count), 0
    )) |>
    dplyr::ungroup()
  return(ss)
}
