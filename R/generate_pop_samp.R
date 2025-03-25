#' Define active infected cell growth curve
#'
#' Get per-generation logistic growth curve for infected active cells.
#' Note that the simulation will only be run until the final sampling time
#' within n_gens. See [define_sampling_scheme()].
#'
#' @param n_gens number of generations (default: 5000)
#' @param carry_cap carrying capacity for number of infected cells (default: 2000)
#' @param n0 starting infected cell population size (default: 10)
#' @param max_growth_rate maximum infected cell population growth rate
#' (default: 0.3)
#'
#' @return tibble with two columns: day and active cell count
#' @export
#'
#' @examples
#' define_growth_curve()
define_growth_curve <- function(n_gens = 5000, n0 = 10, carry_cap = 2000, max_growth_rate = 0.3) {
  check_define_growth_curve_inputs(n_gens, n0, carry_cap, max_growth_rate)
  n <- n0
  sapply(1:n_gens, function(x) {
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
#' Define which days to sample sequences, and how many sequences to sample
#'
#' @param sampling_frequency_active frequency in days at which to record
#'   sequences from active cells (and counts) (default: 365 days)
#' @param max_samp_active maximum number of cells (and thus sequences) to sample
#'   from active cells in a given day (default: 20 sequences)
#' @param sampling_frequency_latent frequency in days at which to record
#'   sequences from latent cells (default: 365 days)
#' @param max_samp_latent maximum number of cells (and thus sequences) to sample
#'   from latent cells in a given day (default: 20 sequences)
#' @param n_days day to end sampling (default: 3650)
#'
#' @return input growth curve tibble with one additional column
#'   (`n_sample_active`) containing the number of sequences from active cells to
#'   samples
#' @export
#'
#' @examples
#' define_sampling_scheme()
define_sampling_scheme <- function(sampling_frequency_active = 365,
                                   max_samp_active = 20,
                                   sampling_frequency_latent = 365,
                                   max_samp_latent = 20,
                                   n_days = 3650) {
  check_define_sampling_scheme_inputs(sampling_frequency_active, max_samp_active, n_days)
  tibble::tibble(
    day = 0:n_days,
    n_sample_active = ifelse(.data$day %% sampling_frequency_active == 0, max_samp_active, 0),
    n_sample_latent = ifelse(.data$day %% sampling_frequency_latent == 0, max_samp_latent, 0)
  ) |>
    dplyr::filter(.data$n_sample_active != 0 | .data$n_sample_latent != 0)
}
