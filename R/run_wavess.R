#' Run wavess
#'
#' Simulate within-host evolution optionally including recombination (default:
#' on), latency (default: on), and fitness costs (default: off). The three
#' fitness costs that can be simulated are conserved sites, fitness relative to
#' a reference sequence, and antibody-based immune fitness costs. Please note
#' that the default arguments were set with the the HIV ENV gp120 gene in mind.
#' If you'd like to simulate something else, you will likely have to modify
#' certain parameters. However, if you are interested in this gene in
#' particular, you can probably use most of the defaults including the founder
#' and reference sequences provided as examples. However, by default no fitness
#' costs are modeled. To model these, you can use the pre-processing functions
#' (see `vignette('prepare_input_data')`) to generate the relevant inputs. Also,
#' the parameters for latent probabilities are assumed to be small, such that it
#' is unlikely that multiple events (activate, die, proliferate) will occur to a
#' single latent cell in a single (active cell) generation. See
#' `vignette('run_wavess')` for more details about the simulator and input
#' arguments.
#'
#' @param pop_samp Tibble with columns generation, active_cell_count,
#'   n_sample_active. Note that the initial active cell population size in the
#'   first generation must be the same as the number of input founder sequences
#'   (because it simply _is_ the input founder sequences). Can be generated
#'   using the [generate_pop_samp()] function.
#' @param founder_seqs Founder sequence(s) as a character string or a vector of
#'   character strings, for example 'ACATG'. The founder sequence(s) may only
#'   contain the characters ACGT, and no gaps are allowed. When modeling
#'   immune fitness, they are expected to be codon-aligned.
#' @param nt_sub_probs Named matrix of nucleotide substitution probabilities.
#'   Rows are from, columns are to. Can be generated using the
#'   [calc_nt_sub_probs()] function.
#' @param conserved_sites Vector of conserved sites. This can be generated using
#'   the [identify_conserved_sites()] function (default: NULL, i.e. no conserved
#'   sites fitness costs)
#' @param conserved_cost Cost of mutation at conserved site (default: 0.99)
#' @param ref_seq Reference sequence as a character string. A consensus
#'   sequence, that can be used as the reference sequence, can be generated
#'   using the function [identify_conserved_sites()] (default: NULL, i.e. no
#'   fitness cost relative to a reference sequence)
#' @param rep_exp Replicative fitness exponent, only relevant when ref_seq is
#'   not NULL (default: 1) # MAKE THIS CLEARER ONCE WE DECIDE ON A FINAL
#'   DEFINITION
#' @param epitope_locations Tibble of epitope locations and maximum fitness
#'   costs with columns epi_start_nt, epi_end_nt, max_fitness_cost. These
#'   epitopes are expected to be in a protein in the correct reading frame, as
#'   the nucleotide sequences are translated to amino acids to calculate the
#'   immune fitness cost. This epitope location tibble can be generated using
#'   the functions [get_epitope_frequencies()] and [sample_epitopes()].
#'   (default: NULL, i.e. no immune fitness costs)
#' @param seroconversion_time Generation at which seroconversion occurs, only
#'   relevant when epitope_locations is not NULL (default: 30).
#' @param prop_for_imm Proportion of all infected cells that must be infected
#'   with a given sequence for that sequence to be recognized by the immune
#'   system, only relevant when epitope_locations is not NULL (default: 0.01).
#' @param gen_full_potency Number of generations it takes for an immune response
#'   to an epitope to reach full potency, only relevant when epitope_locations
#'   is not NULL (default: 90).
#' @param prob_mut Probability of a mutation at one site in one generation
#'   (default: 3.5e-5)
#' @param prob_recomb Probability of a recombination event at a given site in
#'   one generation (default: 1.4e-5)
#' @param prob_act_to_lat Probability that an active cell becomes latent in a
#'   generation (default: 0.001). Set this to 0 if you don't want to model
#'   latent cell dynamics.
#' @param prob_lat_to_act Probability that a latent cell becomes active in a
#'   generation (default: 0.01)
#' @param prob_lat_prolif Probability that a latent cell proliferates in a
#'   generation (default: 0.01)
#' @param prob_lat_die Probability that a latent cell dies in a generation
#'   (default: 0.01)
#' @param seed Optional seed (default: NULL)
#'
#' @return List including: tibble of counts, and alignment of sampled sequences
#' @export
#'
#' @examples
#' \dontrun{
#' run_wavess(
#'   generate_pop_samp(n_gen = 300), "ATCG",
#'   calc_nt_sub_probs(hiv_env_flt_2022)
#' )
#' }
run_wavess <- function(pop_samp,
                       founder_seqs,
                       nt_sub_probs,
                       prob_mut = 3.5e-5,
                       prob_recomb = 1.4e-5,
                       prob_act_to_lat = 0.001,
                       prob_lat_to_act = 0.01,
                       prob_lat_prolif = 0.01,
                       prob_lat_die = 0.01,
                       conserved_sites = NULL,
                       conserved_cost = 0.99,
                       ref_seq = NULL,
                       rep_exp = 1,
                       epitope_locations = NULL,
                       seroconversion_time = 30,
                       prop_for_imm = 0.01,
                       gen_full_potency = 90,
                       seed = NULL) {
  check_run_wavess_inputs(
    pop_samp, founder_seqs, nt_sub_probs,
    prob_mut, prob_recomb,
    conserved_sites, conserved_cost,
    ref_seq, rep_exp,
    epitope_locations, seroconversion_time,
    prop_for_imm, gen_full_potency,
    prob_act_to_lat, prob_lat_to_act,
    prob_lat_prolif, prob_lat_die,
    seed
  )

  # make sure founder sequences are all uppercase
  founder_seqs <- toupper(founder_seqs)

  # initiate python virtual environment
  agents <- use_python_venv()

  latent <- TRUE
  # no latent cells
  if (prob_act_to_lat == 0) {
    latent_nums <- c(0, 0, 0, 0)
    latent <- FALSE
  }

  if (is.null(ref_seq)) {
    ref_seq <- ""
    replicative_fitness <- 0
  } else {
    replicative_fitness <- 1
    ref_seq <- toupper(ref_seq)
    if (!is.null(conserved_sites)) {
      # mask conserved sites so they aren't included in replicative fitness
      # computation
      ref_seq_str <- strsplit(ref_seq, "")[[1]]
      ref_seq_str[conserved_sites] <- "-"
      ref_seq <- paste0(ref_seq_str, collapse = "")
    }
  }
  if (is.null(conserved_sites)) {
    conserved_fitness <- 0
    conserved_sites <- c()
  } else {
    conserved_fitness <- 1
    # change indexing to 0 because underlying functions are in python
    conserved_sites <- conserved_sites - 1
  }
  if (is.null(epitope_locations)) {
    immune_fitness <- 0
  } else {
    immune_fitness <- 1
    # change start indexing to 0 because underlying functions are in python (but
    # don't change end because Python doesn't include the last position in the
    # slice)
    epitope_locations$epi_start_nt <- epitope_locations$epi_start_nt - 1
  }

  conserved_sites <- as.list(conserved_sites)
  if (!is.null(epitope_locations)) {
    epitope_locations <- apply(epitope_locations, 1, function(x) {
      agents$create_epitope(x[1], x[2], x[3])
    })
  }


  # Set seed
  generator <- agents$set_python_seed(seed)

  # Get nucleotide substitution probabilities in right format
  nucleotides_order <- rownames(nt_sub_probs)
  substitution_probabilities <- unname(lapply(
    data.frame(t(nt_sub_probs)),
    function(x) x
  ))

  # Create host environment and initialize infected cells
  host <- agents$create_host_env(
    as.list(founder_seqs),
    ref_seq, conserved_sites, replicative_fitness, rep_exp,
    as.integer(pop_samp$active_cell_count[1])
  )
  # Last sampled generation (don't have to continue simulation after this)
  last_sampled_gen <- max(pop_samp$generation[pop_samp$n_sample_active != 0])

  # Initialize counts and sequences objects
  counts <- tibble::tibble(
    generation = integer(),
    active_cell_count = integer(),
    latent_cell_count = integer(),
    active_turned_latent = integer(),
    latent_turned_active = integer(),
    latent_died = integer(),
    latent_proliferated = integer(),
    number_mutations = integer(),
    number_dual_inf = integer(),
    mean_fitness_active = numeric(),
    mean_conserved_cost_active = numeric(),
    mean_immune_cost_active = numeric(),
    mean_replicative_cost_active = numeric()
  )

  # put founders at top of file
  founders <- strsplit(founder_seqs, "")
  seqs <- ape::as.DNAbin(strsplit(founder_seqs, ""))
  names(seqs) <- paste0("founder", seq_along(founders))

  if (pop_samp$n_sample_active[1] != 0) {
    # num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate
    latent_nums <- c(0, 0, 0, 0)
    # n_mut, number_recombination
    var_nums <- c(0, 0)
    # mean_fitness_active, mean_conserved_cost_active, mean_immune_cost_active,
    # mean_replicative_cost_active
    fitness <- unlist(reticulate::py_to_r(host$summarize_fitness()))
    counts <- record_counts(counts, 0, host, latent_nums, var_nums, fitness)
    seqs <- c(seqs, sample_viral_sequences(
      0, host,
      pop_samp$n_sample_active[1]
    ))
  }

  # Looping through generations until we sample everything we want
  for (t in 1:last_sampled_gen) {
    # Latent reservoir dynamics only get latent cell dynamics if modeling
    # latency (ALSO CHANGE THIS IN MAIN.PY?)
    if (latent) {
      # num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate
      latent_nums <- unlist(reticulate::py_to_r(
        host$get_next_gen_latent(
          prob_act_to_lat, prob_lat_to_act,
          prob_lat_die, prob_lat_prolif, generator
        )
      ))
    }
    # Productively infected cell dynamics
    # n_mut, number_recombination
    var_nums <- unlist(reticulate::py_to_r(host$get_next_gen_active(
      prob_mut, prob_recomb, pop_samp$active_cell_count[t + 1], t,
      seroconversion_time, nucleotides_order, substitution_probabilities,
      conserved_sites, gen_full_potency, conserved_cost, ref_seq,
      prop_for_imm, epitope_locations,
      immune_fitness, conserved_fitness, replicative_fitness, rep_exp, generator
    )))
    # Record events
    if (pop_samp$n_sample_active[t + 1] != 0) {
      fitness <- unlist(reticulate::py_to_r(host$summarize_fitness()))
      counts <- record_counts(counts, t, host, latent_nums, var_nums, fitness)
      seqs <- c(seqs, sample_viral_sequences(
        t, host,
        pop_samp$n_sample_active[t + 1]
      ))
    }
  }

  return(list(counts = counts, seqs = seqs))
}

#' Record counts
#'
#' @param counts Empty tibble or previous counts
#' @param generation Generation
#' @param host Host environment (python)
#' @param latent_nums Numbers related to latency (num_to_make_latent,
#'   num_to_activate, num_to_die, num_to_proliferate)
#' @param var_nums Numbers related to variation (n_mut, number_recombination)
#' @param fitness Numbers related to fitness (mean_fitness_active,
#'   mean_conserved_cost_active, mean_immune_cost_active,
#'   mean_replicative_cost_active)
#'
#' @return Updated counts tibble with additional row added for generation
#' @noRd
record_counts <- function(counts, generation, host, latent_nums, var_nums,
                          fitness) {
  counts |>
    tibble::add_row(
      generation = generation,
      active_cell_count = length(host$C),
      latent_cell_count = length(host$L),
      # num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate
      active_turned_latent = latent_nums[1],
      latent_turned_active = latent_nums[2],
      latent_died = latent_nums[3],
      latent_proliferated = latent_nums[4],
      # n_mut, number_recombination
      number_mutations = var_nums[1],
      number_dual_inf = var_nums[2],
      # mean_fitness_active, mean_conserved_cost_active,
      # mean_immune_cost_active, mean_replicative_cost_active
      mean_fitness_active = fitness[1],
      mean_conserved_cost_active = fitness[2],
      mean_immune_cost_active = fitness[3],
      mean_replicative_cost_active = fitness[4]
    )
}


#' Sample viral sequences
#'
#' @param generation Generation (to include in sequence name)
#' @param host Host environment (python)
#' @param n_to_samp Number of cells to sample
#'
#' @return
#' Sampled sequences in [ape::DNAbin] format.
#' @noRd
sample_viral_sequences <- function(generation, host, n_to_samp) {
  # subtract 1 because python indexes at 0
  sampled_cells <- sample(seq_along(host$C), n_to_samp) - 1
  seqs <- lapply(sampled_cells, function(x) {
    strsplit(reticulate::py_to_r(host$C[[x]]$infecting_virus$nuc_sequence),
      split = ""
    )[[1]]
  })
  names(seqs) <- lapply(sampled_cells, function(x) {
    paste0(
      "gen_", generation, "_cell_", x,
      "_ic_", host$C[[x]]$infecting_virus$immune_fitness_cost,
      "_cc_", host$C[[x]]$infecting_virus$conserved_fitness_cost,
      "_rc_", host$C[[x]]$infecting_virus$replicative_fitness_cost,
      "_f_", host$C[[x]]$infecting_virus$fitness
    )
  })
  return(ape::as.DNAbin(seqs))
}
