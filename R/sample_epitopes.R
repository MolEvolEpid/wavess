#' Sample epitopes
#'
#' Sample epitopes based on epitope probabilities.
#' Note that the positions returned assume that the start of the amino acid
#' sequence is also the start of the founder sequence in the simulation.
#'
#' @param epitope_probabilities Epitope probability tibble as output by
#' `[get_epitope_frequencies()]`, including columns
#' `aa_position` and `epitope_probability`
#' @param start_aa_pos Starting amino acid position to consider for epitopes
#' (default: 1)
#' @param end_aa_pos Ending amino acid position to consider for epitopes
#' (default: NULL, i.e. all positions)
#' @param num_epitopes Number of epitopes to sample
#' @param aa_epitope_length Amino acid epitope length
#' @param max_fit_cost Maximum fitness cost of an epitope,
#' must be between 0 and 1 (default: 0.4)
#' **note that the model output is very sensitive to this parameter**
#' @param cost_type "linear" or "random"; linear returns max fitness costs
#' distributed linearly between 0 and `max_fit_cost` (not including 0, but including
#' `max_fit_cost`), random returns one epitope with `max_fit_cost` and all other
#' epitopes with max fitness costs randomly selected from a uniform distribution
#' between 0 and `max_fit_cost` (default: "linear")
#' @param max_resamples Maximum number of resampling events to attempt;
#' this is to prevent an infinite loop (default: 100)
#' @param ref_founder_map Output from `[map_ref_founder()]`, including reference
#' and founder positions (`ref_pos` and `founder_pos`).
#'
#' @return tibble with the `num_epitopes` rows and the following columns:
#' - `epi_start_nt`: nucleotide epitope start position
#' - `epi_end_nt`: nucleotide epitope end position
#' - `max_fitness_cost`: maximum fitness cost for that epitope
#' @export
#'
#' @examples
#' sample_epitopes(get_epitope_frequencies(env_features$position))
sample_epitopes <- function(epitope_probabilities,
                            start_aa_pos = 1,
                            end_aa_pos = NULL,
                            num_epitopes = 10,
                            aa_epitope_length = 10,
                            max_fit_cost = 0.4,
                            cost_type = 'linear',
                            max_resamples = 100,
                            ref_founder_map = NULL){
  check_sample_epitopes_inputs(epitope_probabilities, start_aa_pos, end_aa_pos,
                               num_epitopes, aa_epitope_length,
                               max_fit_cost, cost_type,
                               max_resamples, ref_founder_map)
  if(is.null(end_aa_pos)){
    end_aa_pos <- max(epitope_probabilities$aa_position)
  }
  # max fitness cost for each epitope
  if(cost_type == 'linear'){
    max_fit_costs <- seq(0, max_fit_cost, length.out = num_epitopes+1)[2:(num_epitopes+1)]
  }else if(cost_type == 'random'){
    max_fit_costs <- c(max_fit_cost, stats::runif(num_epitopes-1, 0, max_fit_cost))
  }
  # Draw num_epitopes positions
  start_pos <- c()  # list of epitope start positions
  all_pos <- c()  # list of all epitope positions
  n_resamples <- 0
  while(length(start_pos) < num_epitopes & n_resamples < max_resamples){
    mid <- sample(epitope_probabilities$aa_position, 1,
                  prob = epitope_probabilities$epitope_probability)
    start <- mid - floor(aa_epitope_length/2)
    if(aa_epitope_length%%2 == 0){ # eptiope has even length
      end <- mid + floor(aa_epitope_length/2) - 1
    }else{
      end <-  mid + floor(aa_epitope_length/2)
    }
    # Don't want overlapping epitopes,
    # and can't start before start of sequence or end after end of sequence
    if(start %in% all_pos | end %in% all_pos | start < start_aa_pos | end > end_aa_pos){
      n_resamples <- n_resamples + 1
    }else{
      start_pos <- c(start_pos, start)
      all_pos <- c(all_pos, start:(end+1))
    }
  }
  if(n_resamples > 0){
    if(n_resamples >= max_resamples){
      stop('Too many resamples required.
            Try increasing `max_resamples`, or decreasing `num_epitopes` or
            `aa_epitope_length`')
    }else{
      message(n_resamples, " resamples required")
    }
  }
  epitopes <- tibble::tibble(epi_start_nt = start_pos*3,
                 epi_end_nt = (start_pos+10)*3,
                 max_fitness_cost = max_fit_costs)
  if(!is.null(ref_founder_map)){
    epitopes <- convert_ref_to_founder_epitopes(epitopes, ref_founder_map)
  }
  return(epitopes)
}

#' Get epitope frequencies
#'
#' @param epitope_positions numeric vector of **amino acid** positions of identified epitopes
#' (e.g. the env features from the [LANL HIV database](https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/download_db.comp))
#'
#' @return tibble with the following columns:
#' - `aa_position`: amino acid position
#' - `n_features`: number of features at that position
#' - `epitope_probability`: estimated probability of an epitope at that position
#' @export
#'
#' @examples
#' get_epitope_frequencies(env_features$position)
get_epitope_frequencies <- function(epitope_positions){
  check_is_pos(epitope_positions)
  # get epitope probability for each site
  tibble::tibble(aa_position = min(epitope_positions):max(epitope_positions)) |>
    dplyr::left_join(tibble::enframe(epitope_positions) |>
                # count up the number of features for each position
                dplyr::group_by(aa_position=epitope_positions) |>
                dplyr::tally(name = 'n_features'), by = 'aa_position') |>
    # make any positions with no features 0
    dplyr::mutate(n_features = ifelse(is.na(.data$n_features),
                                      0, .data$n_features),
           # normalize to get probabilities
           epitope_probability = .data$n_features/sum(.data$n_features))
}

#' Convert reference epitope locations to founder epitope locations
#'
#' @param ref_epitopes Output from `[sample_epitopes()]`, including start and end
#' nt positions of the epitopes (`epi_start_nt` and `epi_end_nt`)
#' @inheritParams sample_epitopes
#'
#' @return Tibble with epitope positions relative to the founder, with the
#' same columns as output by `[sample_epitopes()]`
convert_ref_to_founder_epitopes <- function(ref_epitopes, ref_founder_map){
  # internal function so don't include checks right now...
  ref_epitopes |>
    dplyr::rename(ref_pos_start = 'epi_start_nt',
                  ref_pos_end = 'epi_end_nt') |>
    dplyr::left_join(ref_founder_map |>
                       # if there is a deletion in the HXB2 sequence,
                       # call it the nearest previous position
                       tidyr::fill('founder_pos', .direction = 'down'),
                     by = c('ref_pos_start' = 'ref_pos')) |>
    dplyr::rename(epi_start_nt = 'founder_pos') |>
    dplyr::mutate(epi_end_nt = .data$epi_start_nt + (.data$ref_pos_end - .data$ref_pos_start)) |>
    dplyr::select('epi_start_nt', 'epi_end_nt', 'max_fitness_cost')
}

