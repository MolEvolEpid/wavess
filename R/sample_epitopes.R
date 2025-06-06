#' Sample epitopes
#'
#' Sample epitopes based on epitope probabilities. Note that the positions
#' returned assume that the start of the amino acid sequence is also the start
#' of the founder sequence in the simulation. We also assume that there are no
#' frameshift mutations in the founder sequence.
#'
#' @param epitope_probabilities Epitope probability tibble as output by
#'   [get_epitope_frequencies()], including columns `aa_position` and
#'   `epitope_probability`. `aa_position` should be indexed at 0
#' @param start_aa_pos Starting amino acid position to consider for epitopes,
#'   indexed at 0 (default: 0, i.e. the first position)
#' @param end_aa_pos Ending amino acid position to consider for epitopes,
#'   indexed at 0 (default: NULL, i.e. through the final position in
#'   `epitope_probabilities$aa_position`)
#' @param num_epitopes Number of epitopes to sample
#' @param aa_epitope_length Amino acid epitope length
#' @param max_fit_cost Maximum fitness cost of an epitope, must be in the range
#'   [0,1) where 0 indicates no cost. 1, which indicates no ability to
#'   survive, is not allowed (default: 0.3)
#' **note that the model output is very sensitive to this parameter**
#' @param max_resamples Maximum number of resampling events to attempt; this is
#'   to prevent an infinite loop (default: 100)
#' @param ref_founder_map Output from [map_ref_founder()], including
#'   *nucleotide* reference and founder positions (`ref_pos` and `founder_pos`).
#'   **NOTE:** The reference positions here, if they were converted to amino
#'   acid positions, are expected to match with the reference positions in
#'   `epitope_probabilities`. Further, we assume that the founder indices align
#'   with the founder sequence positions to be used in the simulation (default:
#'   NULL)
#'
#' @return tibble with the `num_epitopes` rows and the following columns:
#' - `epi_start_nt`: nucleotide epitope start position
#' - `epi_end_nt`: nucleotide epitope end position
#' - `max_fitness_cost`: maximum fitness cost for that epitope
#' @export
#'
#' @examples
#' sample_epitopes(get_epitope_frequencies(env_features$Position - 1))
sample_epitopes <- function(epitope_probabilities,
                            start_aa_pos = 0,
                            end_aa_pos = NULL,
                            num_epitopes = 10,
                            aa_epitope_length = 10,
                            max_fit_cost = 0.3,
                            max_resamples = 100,
                            ref_founder_map = NULL) {
  check_sample_epitopes_inputs(
    epitope_probabilities, start_aa_pos, end_aa_pos,
    num_epitopes, aa_epitope_length,
    max_fit_cost,
    max_resamples, ref_founder_map
  )
  if (is.null(end_aa_pos)) {
    end_aa_pos <- max(epitope_probabilities$aa_position)
  }
  # max fitness cost for each epitope
  max_fit_costs <- seq(0, max_fit_cost, length.out = num_epitopes + 1)[2:(num_epitopes + 1)]
  # Draw num_epitopes positions
  start_pos <- c() # list of epitope start positions
  all_pos <- c() # list of all epitope positions
  n_resamples <- 0
  while (length(start_pos) < num_epitopes && n_resamples < max_resamples) {
    mid <- sample(epitope_probabilities$aa_position, 1,
      prob = epitope_probabilities$epitope_probability
    )
    start <- mid - floor(aa_epitope_length / 2)
    if (aa_epitope_length %% 2 == 0) { # eptiope has even length
      end <- mid + floor(aa_epitope_length / 2) - 1
    } else {
      end <- mid + floor(aa_epitope_length / 2)
    }
    # Don't want overlapping epitopes,
    # and can't start before start of sequence or end after end of sequence
    if (start %in% all_pos || end %in% all_pos || start < start_aa_pos ||
      end > end_aa_pos) {
      n_resamples <- n_resamples + 1
    } else {
      start_pos <- c(start_pos, start)
      all_pos <- c(all_pos, start:(end + 1))
    }
  }
  if (n_resamples > 0) {
    if (n_resamples >= max_resamples) {
      stop("Too many resamples required.
            Try increasing `max_resamples`, or decreasing `num_epitopes` or
            `aa_epitope_length`")
    } else {
      message(n_resamples, " resamples required")
    }
  }
  # randomize epitope start positions
  start_pos <- sample(start_pos)
  if (is.null(ref_founder_map)) {
    epitopes <- tibble::tibble(
      # multiply by 3 to get start of amino acid (indexed at 0)
      epi_start_nt = start_pos * 3,
      epi_end_nt = (start_pos + aa_epitope_length) * 3,
      max_fitness_cost = max_fit_costs
    )
  } else {
    epitopes <- reindex_epitopes(
      start_pos, aa_epitope_length, max_fit_costs,
      ref_founder_map
    )
  }
  return(epitopes)
}

#' Get epitope frequencies
#'
#' @param epitope_positions numeric vector of **amino acid** positions of
#'   identified epitopes (e.g. the env features from the [LANL HIV
#'   database](https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/download_db.comp))
#'
#' @return tibble with the following columns:
#' - `aa_position`: amino acid position
#' - `n_features`: number of features at that position
#' - `epitope_probability`: estimated probability of an epitope at that position
#' @export
#'
#' @examples
#' get_epitope_frequencies(env_features$Position - 1) # subtract 1 to index at 0
get_epitope_frequencies <- function(epitope_positions) {
  check_is_pos(epitope_positions, "epitope_positions", ok0 = TRUE)
  # get epitope probability for each site
  tibble::tibble(aa_position = min(epitope_positions):max(epitope_positions)) |>
    dplyr::left_join(tibble::enframe(epitope_positions) |>
      # count up the number of features for each position
      dplyr::group_by(aa_position = epitope_positions) |>
      dplyr::tally(name = "n_features"), by = "aa_position") |>
    # make any positions with no features 0
    dplyr::mutate(
      n_features = ifelse(is.na(.data$n_features),
        0, .data$n_features
      ),
      # normalize to get probabilities
      epitope_probability = .data$n_features / sum(.data$n_features)
    )
}

#' Convert reference epitope locations to founder epitope locations
#'
#' @param start_pos Vector of starting amino acid positions of the epitopes
#' @param max_fit_costs Vector of maximum fitness cost for each epitope
#' @inheritParams sample_epitopes
#'
#' @return Tibble with epitope positions relative to the founder, with the
#' same columns as output by [sample_epitopes()]
#' @noRd
reindex_epitopes <- function(start_pos, aa_epitope_length, max_fit_costs,
                             ref_founder_map) {
  end_pos <- start_pos * 3
  not_in_map <- c(
    start_pos[!(start_pos * 3) %in% ref_founder_map$ref_pos],
    end_pos[!end_pos %in% ref_founder_map$ref_pos]
  )
  if (length(not_in_map)) {
    stop(
      'Not all reference epitope start and end positions are in ",
         "ref_founder_map: ',
      paste0(not_in_map, collapse = ",")
    )
  }
  tibble::tibble(ref_start_pos = start_pos) |>
    dplyr::left_join(
      ref_founder_map |>
        # convert to amino acid positions
        dplyr::mutate(
          ref_start_pos = ceiling(.data$ref_pos / 3),
          founder_start_pos = ceiling(.data$founder_pos / 3)
        ) |>
        # if there is a deletion in the HXB2 sequence,
        # call it the nearest previous position
        tidyr::fill("founder_start_pos", .direction = "down") |>
        dplyr::select("ref_start_pos", "founder_start_pos") |>
        dplyr::distinct() |>
        dplyr::group_by(.data$ref_start_pos) |>
        dplyr::slice_min(.data$founder_start_pos),
      by = "ref_start_pos"
    ) |>
    # multiply by 3 to get start of amino acid (indexed at 0)
    dplyr::mutate(
      epi_start_nt = .data$founder_start_pos * 3,
      epi_end_nt = (.data$founder_start_pos + aa_epitope_length) * 3,
      max_fitness_cost = max_fit_costs
    ) |>
    dplyr::select("epi_start_nt", "epi_end_nt", "max_fitness_cost")
}
