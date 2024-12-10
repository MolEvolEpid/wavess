#' Run wavess
#'
#' Simulate within-host evolution optionally including recombination (default:
#' on), latency (default: on), and fitness costs (default: off). The three
#' fitness costs that can be simulated are conserved sites, fitness relative to
#' a reference sequence, and antibody-based immune fitness costs. Nucleotide
#' positions for conserved and immune fitness are expected to be indexed at 0.
#' Please note that the default arguments were set with the the HIV ENV gp120
#' gene in mind. If you'd like to simulate something else, you will likely have
#' to modify certain parameters. However, if you are interested in this gene in
#' particular, you can probably use most of the defaults including the founder
#' and reference sequences provided as examples. However, by default no fitness
#' costs are modeled. We recommend including fitness to get a realistic model
#' output. To model these, you can use the pre-processing functions (see
#' `vignette('prepare_input_data')`) to generate the relevant inputs. Also, the
#' parameters for latent probabilities are assumed to be small, such that it is
#' unlikely that multiple events (activate, die, proliferate) will occur to a
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
#'   contain the characters ACGT, and no gaps are allowed. When modeling immune
#'   fitness, they are expected to be codon-aligned.
#' @param q Nucleotide substitution rate matrix Q with rows and columns named as
#'   the nucleotides ACGT. Rows are from, columns are to. Can be generated using
#'   the [estimate_q()] function (default: `hiv_q_mat`).
#' @param conserved_sites Vector of conserved bases named by position in the
#'   founder sequence (indexed at 0). This can be generated using the
#'   [identify_conserved_sites()] function (default: NULL, i.e. no conserved
#'   sites fitness costs)
#' @param conserved_cost Cost of mutation at conserved site, must be in the
#'   range [0,1) where 0 indicates no cost. 1, which indicates no ability to
#'   survive, is not allowed (default: 0.99)
#' @param ref_seq Reference sequence as a character string, which denotes the
#'   "most fit" virus from a replicative perspective. A consensus sequence, that
#'   can be used as the reference sequence, can be generated using the function
#'   [identify_conserved_sites()] (default: NULL, i.e. no fitness cost relative
#'   to a reference sequence)
#' @param replicative_cost Replicative fitness cost, only relevant when ref_seq
#'   is not NULL, must be in the range [0,1) where 0 indicates no cost. 1, which
#'   indicates no ability to survive, is not allowed (default: 0.001)
#' @param epitope_locations Tibble of epitope locations and maximum fitness
#'   costs with columns epi_start_nt, epi_end_nt, max_fitness_cost. These
#'   epitopes are expected to be indexed at 0 and in a protein in the correct
#'   reading frame, as the nucleotide sequences are translated to amino acids to
#'   calculate the immune fitness cost. The maximum fitness cost must be in the
#'   range [0,1) where 0 indicates no cost. 1, which indicates no ability to
#'   survive. This epitope location tibble can be generated using the functions
#'   [get_epitope_frequencies()] and [sample_epitopes()]. (default: NULL, i.e.
#'   no immune fitness costs)
#' @param immune_start_time Generation to start checking for an immune response,
#'   only relevant when epitope_locations is not NULL (default: 0, but note that
#'   the immune response will not actually start until there are at least
#'   `n_for_imm` cells in the active population).
#' @param n_for_imm Proportion of all infected cells that must be infected
#'   with a given sequence for that sequence to be recognized by the immune
#'   system, only relevant when epitope_locations is not NULL (default: 100).
#' @param gen_full_potency Number of generations it takes for an immune response
#'   to an epitope to reach full potency, only relevant when epitope_locations
#'   is not NULL (default: 75).
#' @param mut_rate Mutation rate per-site, per-generation (default: 3.5e-5)
#' @param recomb_rate Recombination rate per-site, per-generation (default:
#'   1.4e-5)
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
#' run_wavess(generate_pop_samp(n_gen = 300), rep("ATCG", 10))
#' }
run_wavess <- function(pop_samp,
                       founder_seqs,
                       mut_rate = 3.5e-5,
                       q = wavess::hiv_q_mat,
                       recomb_rate = 1.4e-5,
                       prob_act_to_lat = 0.001,
                       prob_lat_to_act = 0.01,
                       prob_lat_prolif = 0.01,
                       prob_lat_die = 0.01,
                       conserved_sites = NULL,
                       conserved_cost = 0.99,
                       ref_seq = NULL,
                       replicative_cost = 0.001,
                       epitope_locations = NULL,
                       n_for_imm = 100, # should the default be max(pop_samp$active_cell_count)*0.05?
                       gen_full_potency = 75,
                       immune_start_time = 0,
                       seed = NULL) {
  check_run_wavess_inputs(
    pop_samp, founder_seqs, q,
    mut_rate, recomb_rate,
    conserved_sites, conserved_cost,
    ref_seq, replicative_cost,
    epitope_locations, immune_start_time,
    n_for_imm, gen_full_potency,
    prob_act_to_lat, prob_lat_to_act,
    prob_lat_prolif, prob_lat_die,
    seed
  )

  # initiate python virtual environment
  agents <- use_python_venv()

  # convert rates to probabilities
  prob_mut <- rate_to_probability(mut_rate)
  prob_recomb <- rate_to_probability(recomb_rate)

  # make sure founder sequences are all uppercase and in list format
  founder_seqs <- as.list(toupper(founder_seqs))
  # TODO - MAKE USER INPUT NAMED VECTOR
  names(founder_seqs) <- paste0("founder", seq_along(founder_seqs) - 1)

  # convert q matrix to substitution probabilities
  # TODO - CHANGE THIS TO BE PYTHON FUNCTION?
  nt_sub_probs <- calc_nt_sub_probs_from_q(q, mut_rate)
  # Get nucleotide substitution probabilities in right format
  nucleotides_order <- rownames(nt_sub_probs)
  substitution_probabilities <- unname(lapply(
    data.frame(t(nt_sub_probs)),
    function(x) x
  ))

  if (is.null(conserved_sites)) {
    conserved_sites <- reticulate::dict()
  } else {
    conserved_sites <- as.list(toupper(conserved_sites))
  }
  if (is.null(ref_seq)) {
    ref_seq <- ""
  } else {
    ref_seq <- toupper(ref_seq)
    if (!is.null(conserved_sites)) {
      prep_out <- reticulate::py_to_r(agents$prep_ref_conserved(founder_seqs, ref_seq, conserved_sites))
      ref_seq <- prep_out[[1]]
      conserved_sites <- prep_out[[2]]
    }
  }
  if (!is.null(epitope_locations)) {
    epitope_locations <- apply(epitope_locations, 1, function(x) {
      agents$create_epitope(x[1], x[2], x[3])
    })
  }

  # Set seed
  generator <- agents$set_python_seed(seed)

  # Create host environment and initialize infected cells
  host <- agents$create_host_env(
    founder_seqs,
    ref_seq, replicative_cost,
    as.integer(pop_samp$active_cell_count[1])
  )
  # Last sampled generation (don't have to continue simulation after this)
  last_sampled_gen <- max(pop_samp$generation[pop_samp$n_sample_active != 0])

  # Simulate within-host evolution
  out <- reticulate::py_to_r(host$loop_through_generations(
    pop_samp$active_cell_count,
    pop_samp$n_sample_active,
    last_sampled_gen,
    founder_seqs, nucleotides_order, substitution_probabilities,
    prob_mut, prob_recomb,
    prob_act_to_lat, prob_lat_to_act, prob_lat_die, prob_lat_prolif,
    conserved_sites, conserved_cost, ref_seq, replicative_cost,
    epitope_locations, immune_start_time, n_for_imm, gen_full_potency,
    generator
  ))

  # Fix up output
  names(out) <- c("counts", "fitness", "seqs")
  out$counts <- out$counts |> dplyr::bind_rows()
  out$fitness <- out$fitness |> dplyr::bind_rows()
  out$seqs <- ape::as.DNAbin(t(sapply(out$seqs, function(x) strsplit(x, split = "")[[1]])))

  # Return output
  return(out)
}
