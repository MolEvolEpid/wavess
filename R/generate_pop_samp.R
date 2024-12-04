#' Generate population growth curve and sampling scheme
#'
#' @param n_gen final sampling generation (default: 3000)
#' @param carry_cap carrying capacity (default: 2000)
#' @param n0 starting population size (default: 10)
#' @param max_growth_rate maximum infected cell population growth rate
#' (default: 0.3)
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
                              n0 = 10,
                              max_growth_rate = 0.3,
                              sampling_frequency = 300,
                              max_samp = 20) {
  check_generate_pop_samp_inputs(
    n_gen, carry_cap, n0, max_growth_rate,
    sampling_frequency, max_samp
  )
  define_growth_curve(n_gen, n0, carry_cap, max_growth_rate) |>
    define_sampling_scheme(sampling_frequency, max_samp)
}

# Get logistic growth curve for infected active cells

#' Define growth curve
#'
#' @inheritParams generate_pop_samp
#'
#' @return tibble with two columns: generation and active cell count
#' @noRd
define_growth_curve <- function(n_gen = 3000, n0 = 1, carry_cap = 2000, max_growth_rate = 1) {
  n <- n0
  sapply(1:n_gen, function(x) {
    n <<- min(n + max_growth_rate * n * (carry_cap - n) / carry_cap, carry_cap)
    return(c(x, n))
  }) |>
    t() |>
    tibble::as_tibble(.name_repair = "unique") |>
    suppressMessages() |>
    dplyr::rename(generation = "...1", active_cell_count = "...2") |>
    tibble::add_row(generation = 0, active_cell_count = n0, .before = 1) |>
    dplyr::mutate(active_cell_count = ceiling(.data$active_cell_count))
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
