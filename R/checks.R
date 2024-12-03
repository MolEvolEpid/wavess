#' Check generate_pop_samp
#'
#' @inheritParams generate_pop_samp
#'
#' @return error if inputs are incorrect
#' @noRd
check_generate_pop_samp_inputs <- function(n_gen,
                                           carry_cap,
                                           n0,
                                           max_growth_rate,
                                           sampling_frequency, max_samp) {
  check_is_pos(n_gen, "n_gen")
  check_is_pos(carry_cap, "carry_cap")
  check_is_numeric(n0, "n0")
  if (n0 > carry_cap) {
    stop("n0 must be a number \u2264carry_cap, but is ", n0)
  }
  check_is_pos(max_growth_rate, "max_growth_rate")
  check_is_pos(sampling_frequency, "sampling_frequency")
  check_is_pos(max_samp, "max_samp")
  if (sampling_frequency > n_gen) {
    stop(
      "sampling_frequency must be a number \u2264maximum generation, but is ",
      sampling_frequency
    )
  }
}

#' Check get_seq_pos inputs
#'
#' @inheritParams get_seq_pos
#'
#' @return error if inputs are incorrect
#' @noRd
check_get_seq_pos_inputs <- function(aln_df, col_name) {
  check_is_df(aln_df, "aln_df")
  check_is_string(col_name, "col_name")
}

#' Check map_ref_founder inputs
#'
#' @inheritParams map_ref_founder
#'
#' @return error if inputs are incorrect
#' @noRd
check_map_ref_founder_inputs <- function(aln, ref, founder) {
  check_is_dnabin(aln, "aln")
  check_is_string(ref, "ref")
  check_is_string(founder, "founder")
  check_name_in_alignment(aln, ref, "aln", "ref")
  check_name_in_alignment(aln, founder, "aln", "founder")
}

#' Check find_consensus_inputs
#'
#' @inheritParams find_consensus
#'
#' @return error if inputs are incorrect
#' @noRd
check_find_consensus_inputs <- function(aln, founder, ref, founder_aln) {
  check_is_dnabin(aln, "aln")
  check_is_string(founder, "founder")
  if (is.null(ref) && is.null(founder_aln)) {
    check_name_in_alignment(aln, founder, "aln", "founder")
  } else if (!is.null(ref)) {
    check_is_string(ref, "ref")
    if (is.null(founder_aln)) {
      stop("When `ref` is specified, `founder_aln` must also be specified")
    }
  } else if (!is.null(founder_aln)) {
    check_is_dnabin(founder_aln, "founder_aln")
    if (is.null(ref)) {
      stop("When `founder_aln` is specified, `ref` must also be specified")
    }
  }
  if (!is.null(ref) && !is.null(founder_aln)) {
    check_name_in_alignment(aln, ref, "aln", "ref")
    check_name_in_alignment(founder_aln, ref, "founder_aln", "ref")
    check_name_in_alignment(founder_aln, founder, "founder_aln", "founder")
  }
}

#' Check identify_conserved_sites inputs
#'
#' @inheritParams identify_conserved_sites
#'
#' @return error if inputs are incorrect
#' @noRd
check_conserved_sites_inputs <- function(aln, founder, thresh, ref,
                                         founder_aln) {
  check_find_consensus_inputs(aln, founder, ref, founder_aln)
  check_is_numeric(thresh)
  if (thresh < 0 || thresh > 1) {
    stop("`thresh` must be a number between 0 and 1 inclusive")
  }
}


#' Check if name is in alignment
#'
#' @param aln Alignment in DNAbin format
#' @param name Name of sequence to check
#' @param aln_name Name of alignment
#' @param name_name Name of name
#'
#' @return error if name isn't in the alignment
#' @noRd
check_name_in_alignment <- function(aln, name, aln_name, name_name) {
  check_is_dnabin(aln, aln_name)
  check_is_string(name, name_name)
  # added this because of a weird issue with devtools::check()...
  aln <- ape::as.matrix.DNAbin(aln)
  if (!name %in% labels(aln)) {
    stop(
      name_name, " must be the name of a sequence in ", aln_name, ", but is ",
      name
    )
  }
}

#' Check if a variable is numeric
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not numeric
#' @noRd
check_is_numeric <- function(x, var_name) {
  if (!is.numeric(x)) {
    stop(var_name, " must be numeric, but is a ", class(x))
  }
}

#' Check if a variable is positive
#'
#' @param ok0 whether 0 is okay (default: FALSE)
#'
#' @inheritParams check_is_numeric
#'
#' @return error if variable is not positive
#' @noRd
check_is_pos <- function(x, var_name, ok0 = FALSE) {
  check_is_numeric(x, var_name)
  if (ok0) {
    if (any(x < 0)) {
      stop(var_name, " must be a number >= 0, but is ", x)
    }
  } else {
    if (any(x <= 0)) {
      stop(var_name, " must be a positive number(s), but is ", x)
    }
  }
}

#' Check if a variable is a dataframe or tibble
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not a dataframe
#' @noRd
check_is_df <- function(x, var_name) {
  if (!is.data.frame(x)) {
    stop(var_name, " must be a data frame or tibble, but is a ", class(x))
  }
}

#' Check if a variable is a matrix
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not a matrix
#' @noRd
check_is_matrix <- function(x, var_name) {
  if (!is.matrix(x)) {
    stop(var_name, " must be a matrix, but is a ", class(x))
  }
}

#' Check if a variable is a string
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not a string
#' @noRd
check_is_string <- function(x, var_name) {
  if (!is.character(x)) {
    stop(var_name, " must be a string, but is a ", class(x))
  }
}

#' Check if a variable is between 0 and 1 inclusive
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not between 0 and 1 inclusive
#' @noRd
check_is_0to1 <- function(x, var_name, ok1 = TRUE) {
  check_is_numeric(x, var_name)
  check_bool <- x < 0 || x > 1
  if (!ok1) {
    check_bool <- x < 0 || x >= 1
  }
  if (check_bool) {
    stop(var_name, " must be in the range [0,1), but is ", x)
  }
}

#' Check if a variable is a DNAbin object
#'
#' @param x variable to check
#' @param var_name name of variable
#'
#' @return error if variable is not a DNAbin object
#' @noRd
check_is_dnabin <- function(x, var_name) {
  class_x <- class(x)
  if (!class_x == "DNAbin") {
    stop(
      var_name,
      " must be of the `ape` class `DNAbin` (e.g. an alignment read in using ",
      "`ape::read.FASTA()` or `ape::read.dna()` but is ",
      class_x
    )
  }
}


#' Check if a variable is a phylo object
#'
#' @param x variable to check
#' @param var_name name of variable
#' @param rooted whether the tree must be rooted
#'
#' @return error if variable is not a DNAbin object
#' @noRd
check_is_phylo <- function(x, var_name, rooted = FALSE) {
  class_x <- class(x)
  if (!class_x == "phylo") {
    stop(
      var_name,
      " must be of the `ape` class `phylo` (e.g. an tree read in using ",
      "`ape::read.tree()` but is ",
      class_x
    )
  }
  if (rooted) {
    if (!ape::is.rooted(x)) {
      stop(var_name, " must be rooted")
    }
  }
}

check_sample_epitopes_inputs <- function(epitope_probabilities,
                                         start_aa_pos, end_aa_pos,
                                         num_epitopes, aa_epitope_length,
                                         max_fit_cost,
                                         max_resamples, ref_founder_map) {
  check_is_df(epitope_probabilities)
  if (!all(c("aa_position", "epitope_probability") %in%
    colnames(epitope_probabilities))) {
    stop("The following columns must be included in `epitope_probabilities`:
         `aa_position`, `epitope_probability`. See `sample_epitopes()`
         documentation for more details")
  }
  check_is_pos(start_aa_pos, "start_aa_pos", ok0 = TRUE)
  if (!is.null(end_aa_pos)) {
    check_is_pos(end_aa_pos, "end_aa_pos")
  }
  check_is_pos(num_epitopes, "num_epitopes")
  check_is_pos(aa_epitope_length, "aa_epitope_length")
  check_is_0to1(max_fit_cost, "max_fit_cost", ok1 = FALSE)
  check_is_pos(max_resamples, "max_resamples")
  if (!is.null(ref_founder_map)) {
    check_is_df(ref_founder_map)
    if (!all(c("ref_pos", "founder_pos") %in% colnames(ref_founder_map))) {
      stop("The following columns must be included in `ref_founder_map`:
         `ref_pos`, `founder_pos`. See `sample_epitopes()`
         documentation for more details")
    } else {
      check_is_numeric(ref_founder_map$ref_pos, "ref_pos in ref_founder_map")
      check_is_numeric(
        ref_founder_map$founder_pos,
        "founder_pos in ref_founder_map"
      )
    }
  }
}

#' Check estimate_q inputs
#'
#' @inheritParams estimate_q
#'
#' @return error if inputs are incorrect
#' @noRd
check_estimate_q_inputs <- function(aln, tr, model, rearrangement) {
  check_is_dnabin(aln, "aln")
  if (!is.null(tr)) {
    check_is_phylo(tr, "tr")
  }
  check_is_string(model, "model")
  check_is_string(rearrangement, "rearrangement")
  if (!rearrangement %in% c("none", "NNI", "stochastic", "ratchet")) {
    stop(
      'Rearrangement must be one of "none", "NNI", "stochastic" or "ratchet", ",
      "but is ',
      rearrangement
    )
  }
}

#' Check check_seq sequence input
#'
#' @param seq character sequence
#' @param chars characters
#' @param seq_type type of sequence
#'
#' @return error if wrong characters in sequence
#' @noRd
check_seq <- function(seq, chars, seq_type) {
  if (!all(strsplit(seq, "")[[1]] %in% chars)) {
    stop(
      seq_type, " must only contain the characters ",
      paste0(chars, collapse = "")
    )
  }
}

#' Check run_wavess inputs
#'
#' @inheritParams run_wavess
#'
#' @return error if wrong inputs
#' @noRd
check_run_wavess_inputs <- function(pop_samp, founder_seqs, q,
                                    mut_rate, recomb_rate,
                                    conserved_sites, conserved_cost,
                                    ref_seq, replicative_cost,
                                    epitope_locations, seroconversion_time,
                                    prop_for_imm, gen_full_potency,
                                    prob_act_to_lat, prob_lat_to_act,
                                    prob_lat_prolif, prob_lat_die,
                                    seed) {
  check_is_df(pop_samp, "pop_samp")
  if (!all(c("generation", "active_cell_count", "n_sample_active") %in%
    colnames(pop_samp))) {
    stop(
      "pop_samp must contain the columns ",
      "generation, active_cell_count, n_sample_active"
    )
  }
  check_is_pos(pop_samp$generation, "pop_samp$generation", TRUE)
  check_is_pos(pop_samp$active_cell_count, "pop_samp$active_cell_count", TRUE)
  if (all(pop_samp$n_sample_active == 0)) {
    stop("you must sample at least one generation")
  }
  check_is_pos(pop_samp$n_sample_active, "pop_samp$n_sample_active", TRUE)
  if (!all(pop_samp$generation == seq_len(nrow(pop_samp)) - 1)) {
    stop(
      "pop_samp$generation must be consecutive numbers from 0 to ",
      "the number of rows in the data"
    )
  }
  if (!all(pop_samp$active_cell_count >= pop_samp$n_sample_active)) {
    stop(
      "pop_samp$active_cell_count must always be ",
      ">= pop_samp$n_sample_active"
    )
  }
  check_is_string(founder_seqs, "founder_seqs")
  if (pop_samp$active_cell_count[1] != length(founder_seqs)) {
    stop("Initial population size must equal the number of founder sequences")
  }
  lapply(as.list(founder_seqs), function(x) {
    check_seq(toupper(x), c("A", "C", "G", "T"), "founder_seqs")
  })
  if (length(unique(sapply(as.list(founder_seqs), nchar))) != 1) {
    stop("All founder sequences must be the same length")
  }
  check_is_matrix(q, "q")
  if (!all(rownames(q) %in% c("A", "C", "G", "T")) ||
    !all(c("A", "C", "G", "T") %in% rownames(q))) {
    stop("q must have rownames A,C,G,T")
  }
  if (!all(colnames(q) %in% c("A", "C", "G", "T")) ||
    !all(c("A", "C", "G", "T") %in% colnames(q))) {
    stop("q must have colnames A,C,G,T")
  }
  if (!all(rownames(q) == colnames(q))) {
    stop("the row names and column names of q must be in the same order")
  }
  sapply(q, function(x) {
    lapply(x, function(y) {
      check_is_numeric(y, "q")
    })
  })
  if (!is.null(conserved_sites)) {
    check_seq(paste0(toupper(conserved_sites), collapse = ""), c("A", "C", "G", "T"), "conserved_sites")
    if (is.null(names(conserved_sites))) {
      stop("conserved_sites must be a named vector")
    }
    check_is_pos(as.numeric(names(conserved_sites)), "names of conserved_sites", ok0 = TRUE)
    if (max(as.numeric(names(conserved_sites))) > nchar(as.list(founder_seqs)[[1]])) {
      stop(
        "the maximum value of conserved_sites is greater than ",
        "the length of the founder sequence"
      )
    }
    check_is_0to1(conserved_cost, "conserved_cost", ok1 = FALSE)
  }
  if (!is.null(ref_seq)) {
    check_is_string(ref_seq, "ref_seq")
    check_seq(toupper(ref_seq), c("A", "C", "G", "T", "-"), "ref_seq")
    if (nchar(as.list(founder_seqs)[1]) != nchar(ref_seq)) {
      stop("ref_seq must be the same length as the founder sequence(s)")
    }
    check_is_0to1(replicative_cost, "replicative_cost", ok1 = FALSE)
  }
  if (!is.null(epitope_locations)) {
    check_is_df(epitope_locations, "epitope_locations")
    if (!all(c("epi_start_nt", "epi_end_nt", "max_fitness_cost") %in%
      colnames(epitope_locations))) {
      stop(
        "epitope_locations must contain the columns ",
        "epi_start_nt, epi_end_nt, max_fitness_cost"
      )
    }
    lapply(epitope_locations$epi_start_nt, function(x) {
      check_is_pos(x, "epitope_locations$epi_start_nt", TRUE)
    })
    lapply(epitope_locations$epi_end_nt, function(x) {
      check_is_pos(x, "epitope_locations$epi_end_nt", FALSE)
    })
    lapply(epitope_locations$max_fitness_cost, function(x) {
      check_is_0to1(x, "epitope_locations$max_fitness_cost", ok1 = FALSE)
    })
    check_is_pos(seroconversion_time, "seroconversion_time", TRUE)
    check_is_0to1(prop_for_imm, "prop_for_imm")
    check_is_pos(gen_full_potency, "gen_full_potency", TRUE)
  }
  check_is_numeric(mut_rate, "mut_rate")
  check_is_numeric(recomb_rate, "recomb_rate")
  check_is_0to1(prob_act_to_lat, "prob_act_to_lat")
  check_is_0to1(prob_lat_to_act, "prob_lat_to_act")
  check_is_0to1(prob_lat_prolif, "prob_lat_prolif")
  check_is_0to1(prob_lat_die, "prob_lat_die")
  if (!is.null(seed)) {
    check_is_numeric(seed, "seed")
  }
}

#' Check extract founder inputs
#'
#' @inheritParams extract_seqs
#'
#' @return error if wrong inputs
#' @noRd
check_extract_seqs_inputs <- function(aln, founder_name, ref_name, start, end) {
  check_name_in_alignment(aln, founder_name, "aln", "founder_name")
  if (!is.null(ref_name)) {
    check_name_in_alignment(aln, ref_name, "aln", "ref_name")
  }
  aln <- as.matrix(aln)
  check_is_pos(start, "start")
  if (!is.null(end)) {
    check_is_pos(end, "end")
    if (end > ncol(aln)) {
      stop("end must be <= the length of the alignment")
    }
  }
}

#' Check slice_aln inputs
#'
#' @inheritParams slice_aln
#'
#' @return error if wrong inputs
#' @noRd
check_slice_aln_inputs <- function(aln, start, end, seqs) {
  lapply(seqs, function(x) {
    check_name_in_alignment(aln, x, "aln", "seqs")
  })
  aln <- as.matrix(aln)
  check_is_pos(start, "start")
  if (!is.null(end)) {
    check_is_pos(end, "end")
    if (end > ncol(aln)) {
      stop("end must be <= the length of the alignment")
    }
  }
}
