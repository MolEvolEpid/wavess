#' Generate population growth curve and sampling scheme
#'
#' @param n_gen final sampling generation (default: 3000)
#' @param carry_cap carrying capacity (default: 2000)
#' @param n0 starting population size (default: 1; used for logistic)
#' @param g50 generation at which the population size reaches 50% of carry_cap
#' (default: 25; used for logistic)
#' @param sampling_frequency frequency in generations at which to record
#' sequences (and counts) (default: 300 generations)
#' @param max_samp maximum number of cells (and thus sequences) to sample in a
#' given generation (default: 20 sequences)
#'
#' @return tibble with two columns:
#' - `generation`: Each generation to be simulated. Here, we consider generation 0 to
#' be the first generation. In the simulation, this generation will consists of
#' the input founder sequences.
#' - `active_cell_count`: Number of active cells in each generation
#' - `n_sample_active` Number of sequences to sample in each generation
#' @export
#'
#' @examples
#' generate_pop_samp()
generate_pop_samp <- function(n_gen = 3000,
                              carry_cap = 2000,
                              n0 = 1,
                              g50 = 25,
                              sampling_frequency = 300,
                              max_samp = 20) {
  check_generate_pop_samp_inputs(
    n_gen, carry_cap, n0, g50,
    sampling_frequency, max_samp
  )
  define_growth_curve(n_gen, carry_cap, n0, g50, 0.5) |>
    define_sampling_scheme(sampling_frequency, max_samp)
}

# Get logistic growth curve for infected active cells

#' Define growth curve
#'
#' @param gen_prop_carry_cap generation at which the population size reaches
#'   prop_carry_cap of carry_cap (default: 25; used for logistic)
#' @param prop_carry_cap proportion of carrying capacity population size at
#'   gen_prop_carry_cap (default: 0.5; used for logistic)
#' @inheritParams generate_pop_samp
#'
#' @return tibble with two columns: generation and active cell count
#' @noRd
define_growth_curve <- function(n_gen = 3000,
                                carry_cap = 2000,
                                n0 = 1,
                                gen_prop_carry_cap = 25,
                                prop_carry_cap = 0.5) {
  gen_df <- tibble::tibble(
    generation = 0:n_gen,
    # infected cell population size over time
    active_cell_count = get_logistic_curve(
      n0, carry_cap, gen_prop_carry_cap,
      prop_carry_cap, n_gen
    )
  )
  return(gen_df)
}

#' Get logistic growth curve
#'
#' Note: to get the equations used in this function, I solved for the growth
#' rate and midpoint based on n0, carry_cap, gen_prop_carry_cap, and
#' prop_carry_cap.
#'
#' @param gen_prop_carry_cap generation at which the population size reaches pK
#'   of carry_cap (default: 25; used for logistic)
#' @param prop_carry_cap proportion of carrying capacity population size at
#'   gen_prop_carry_cap (default: 0.5; used for logistic)
#' @inheritParams generate_pop_samp
#'
#' @return Population size at each generation
#' @noRd
get_logistic_curve <- function(n0, carry_cap, gen_prop_carry_cap,
                               prop_carry_cap, n_gen) {
  # number of cells at seroconversion generation
  n_s <- carry_cap * prop_carry_cap
  # useful constants
  c0 <- log(carry_cap / n0 - 1)
  c_s <- log(carry_cap / n_s - 1)
  # get growth rate
  k <- -(c_s - c0) / gen_prop_carry_cap
  # get midpoint
  x_m <- c0 / k
  return(floor(sapply(0:n_gen, function(x) logistic_fn(x, carry_cap, k, x_m))))
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
logistic_fn <- function(x, carry_cap, growth_rate, midpoint) {
  carry_cap / (1 + exp(-growth_rate * (x - midpoint)))
}

#' Define sampling scheme
#'
#' @param growth_curve output from [define_growth_curve()], or customized growth
#'   curve with columns named `generation` and `active_cell_count`
#' @inheritParams generate_pop_samp
#'
#' @return input growth curve tibble with one additional column
#'   (`n_sample_active`) containing the number of sequences from active cells to
#'   samples
#' @noRd
define_sampling_scheme <- function(growth_curve,
                                   sampling_frequency = 300,
                                   max_samp = 20) {
  ss <- growth_curve |>
    dplyr::rowwise() |>
    # sample up to max_samp sequences every sampling_frequency generations
    dplyr::mutate(n_sample_active = ifelse(.data$generation %%
      sampling_frequency == 0,
    min(max_samp, .data$active_cell_count), 0
    )) |>
    dplyr::ungroup()
  return(ss)
}
