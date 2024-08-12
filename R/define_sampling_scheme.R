#' Define sampling scheme
#'
#' @param growth_curve output from `define_growth_curve()`, or customized
#' growth curve with columns named `generation` and `active_cell_count`
#' @param sampling_frequency frequency in generations at which to record
#' sequences (and counts) (default: 300 generations)
#' @param max_samp maximum number of cells (and thus sequences) to sample in a
#' given generation (default: 20 sequences)
#'
#' @return input growth curve tibble with one additional column (`n_sample_active`)
#' containing the number of sequences from active cells to samples
#' @export
#'
#' @examples
#' define_growth_curve() |> define_sampling_scheme()
define_sampling_scheme <- function(growth_curve,
                                   sampling_frequency = 300,
                                   max_samp = 20){
  check_define_sampling_scheme_inputs(growth_curve, sampling_frequency, max_samp)
  ss <- growth_curve |>
    dplyr::rowwise() |>
    # sample up to max_samp sequences every sampling_frequency generations
    dplyr::mutate(n_sample_active = ifelse(.data$generation%%sampling_frequency == 0,
                                    min(max_samp, .data$active_cell_count), 0)) |>
    dplyr::ungroup()
    return(ss)
}
