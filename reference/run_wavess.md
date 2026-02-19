# Run wavess

Simulate within-host evolution optionally including recombination
(default: on), latency (default: on), and fitness costs (default: off).
The four fitness costs that can be simulated are conserved sites,
fitness relative to a reference sequence, B-cell immune fitness costs,
and T-cell immune fitness costs. Nucleotide positions for conserved and
immune fitness are expected to be indexed at 0. Please note that the
default arguments were set with the the HIV *env* gp120 gene in mind. If
you'd like to simulate something else, you will likely have to modify
certain parameters. However, if you are interested in this gene in
particular, you can probably use most of the defaults including the
founder and reference sequences provided as examples. However, by
default no fitness costs are modeled. We recommend including fitness to
obtain a realistic model output. To model these, you can use the
pre-processing functions (see
[`vignette('prepare_input_data')`](https://molevolepid.github.io/wavess/articles/prepare_input_data.md))
to generate the relevant inputs. Also, the parameters for latent
probabilities are assumed to be small, such that it is unlikely that
multiple events (activate, die, proliferate) will occur to a single
latent cell in a single (active cell) generation. See
[`vignette('run_wavess')`](https://molevolepid.github.io/wavess/articles/run_wavess.md)
for more details about the simulator and input arguments.

## Usage

``` r
run_wavess(
  inf_pop_size,
  samp_scheme,
  founder_seqs,
  generation_time = 1,
  mut_rate = 3e-05,
  q = wavess::calc_q_from_rates(wavess::hiv_mut_rates, mut_rate, generation_time),
  recomb_rate = 1.5e-05,
  act_to_lat = 0.001,
  lat_to_act = 0.01,
  lat_prolif = 0.01,
  lat_die = 0.01,
  conserved_sites = NULL,
  conserved_cost = 0.99,
  ref_seq = NULL,
  replicative_cost = 0.001,
  b_epitope_locations = NULL,
  b_immune_start_day = 0,
  b_n_for_imm = 100,
  b_days_full_potency = 90,
  epitope_locations = NULL,
  n_for_imm = 100,
  days_full_potency = 90,
  immune_start_day = 0,
  t_epitope_locations = NULL,
  t_max_immune_cost = 0.5,
  seed = NULL
)
```

## Arguments

- inf_pop_size:

  Tibble with columns *generation* (starting at day 0) and
  active_cell_count. Note that the initial active cell population size
  on day 0 must be the same as the number of input founder sequences
  (because it simply *is* the input founder sequences). Can be generated
  using the
  [`define_growth_curve()`](https://molevolepid.github.io/wavess/reference/define_growth_curve.md)
  function.

- samp_scheme:

  Tibble with columns *day* and n_sample_active. Rows only need to
  contain the days on which sampling occurrs. Can be generated using the
  [`define_sampling_scheme()`](https://molevolepid.github.io/wavess/reference/define_sampling_scheme.md)
  function.

- founder_seqs:

  Founder sequence(s) as a character string or a vector of character
  strings, for example 'ACATG'. The founder sequence(s) may only contain
  the characters ACGT, and no gaps are allowed. When modeling immune
  fitness, they are expected to be codon-aligned.

- generation_time:

  Amount of time in days it takes a virus to complete one full life
  cycle, from infecting one cell to exiting the cell and infecting the
  next one (default: 1 day). Any inputs that are in days will be
  converted to generations using this number.

- mut_rate:

  Mutation rate per-site, per-*generation* (default: 3.0e-5)

- q:

  Nucleotide substitution rate matrix Q with rows and columns named as
  the nucleotides ACGT. Rows are from, columns are to. Can be generated
  using the
  [`estimate_q()`](https://molevolepid.github.io/wavess/reference/estimate_q.md)
  function. The default is to calculate the Q matrix using estimates of
  per-day rates from nearly neutral sites:
  `wavess::calc_q_from_rates(wavess::hiv_mut_rates,mut_rate,generation_time)`.

- recomb_rate:

  Recombination rate per-breakpoint, per-*generation* (default: 1.5e-5).
  This can be a single number or a numeric vector where each element in
  the vector is a breakpoint-specific recombination rate. If the input
  is a vector, the length of the vector should be one fewer than the
  number of basepairs in the founder sequence. Note that this rate is
  modeled as an effective recombination rate that includes the rate of
  co-infection followed by recombination.

- act_to_lat:

  Per-*day* rate that an active cell becomes latent (default: 0.001).
  Set this to 0 if you don't want to model latent cell dynamics.

- lat_to_act:

  Per-*day* rate that a latent cell becomes active (default: 0.01)

- lat_prolif:

  Per-*day* rate that a latent cell proliferates (default: 0.01)

- lat_die:

  Per-*day* rate that a latent cell dies (default: 0.01)

- conserved_sites:

  Vector of conserved bases named by position in the founder sequence
  (indexed at 0). This can be generated using the
  [`identify_conserved_sites()`](https://molevolepid.github.io/wavess/reference/identify_conserved_sites.md)
  function (default: NULL, i.e. no conserved sites fitness costs)

- conserved_cost:

  Cost of mutation at conserved site, must be in the range \[0,1) where
  0 indicates no cost. 1, which indicates no ability to survive, is not
  allowed (default: 0.99)

- ref_seq:

  Reference sequence as a character string, which denotes the "most fit"
  virus from a replicative perspective. A consensus sequence, that can
  be used as the reference sequence, can be generated using the function
  [`identify_conserved_sites()`](https://molevolepid.github.io/wavess/reference/identify_conserved_sites.md)
  (default: NULL, i.e. no fitness cost relative to a reference sequence)

- replicative_cost:

  Replicative fitness cost, only relevant when ref_seq is not NULL, must
  be in the range \[0,1) where 0 indicates no cost. 1, which indicates
  no ability to survive, is not allowed (default: 0.001)

- b_epitope_locations:

  Tibble of B-cell epitope locations and maximum fitness costs with
  columns epi_start_nt, epi_end_nt, max_fitness_cost. These epitopes are
  expected to be indexed at 0 and in a protein in the correct reading
  frame, as the nucleotide sequences are translated to amino acids to
  calculate the B-cell immune fitness cost. The maximum fitness cost
  must be in the range \[0,1) where 0 indicates no cost. 1, which
  indicates no ability to survive. This epitope location tibble can be
  generated using the functions
  [`get_epitope_frequencies()`](https://molevolepid.github.io/wavess/reference/get_epitope_frequencies.md)
  and
  [`sample_epitopes()`](https://molevolepid.github.io/wavess/reference/sample_epitopes.md).
  (default: NULL, i.e. no B-cell immune fitness costs)

- b_immune_start_day:

  *Day* to start checking for a B-cell immune response, only relevant
  when `b_epitope_locations` is not NULL (default: 0, but note that the
  immune response will not actually start until there are at least
  `b_n_for_imm` cells in the active population).

- b_n_for_imm:

  Number of infected cells that must contain a given sequence for that
  sequence to be recognized by the B-cell immune system, only relevant
  when `b_epitope_locations` is not NULL (default: 100).

- b_days_full_potency:

  Number of *days* it takes for a B-cell immune response to an epitope
  to reach full potency, only relevant when `b_epitope_locations` is not
  NULL (default: 90).

- epitope_locations:

  Deprecated; use `b_epitope_locations` instead. If both
  `epitope_locations` and `b_epitope_locations` are supplied,
  `b_epitope_locations` will be used.

- n_for_imm:

  Deprecated; use `b_n_for_imm` instead.

- days_full_potency:

  Deprecated; use `b_days_full_potency` instead.

- immune_start_day:

  Deprecated; use `b_immune_start_day` instead.

- t_epitope_locations:

  Optional tibble of T-cell epitope locations and escape information
  with columns start (nucleotide start position, indexed at 0),
  days_to_full_potency (days to reach full immune potency for that
  epitope), escape_position (amino acid position within the epitope,
  indexed starting at 1), and recognized_aa (amino acid considered
  recognized by the immune system at that escape position). When
  provided, these are used to compute an additional T-cell immune
  fitness cost. (default: NULL, i.e. no T-cell immune fitness costs)

- t_max_immune_cost:

  Maximum fitness cost per recognized T-cell epitope, must be in the
  range \[0,1) where 0 indicates no cost. 1, which indicates no ability
  to survive, is not allowed (default: 0.5).

- seed:

  Optional seed (default: NULL)

## Value

List including: tibble of counts and mean fitness values, an alignment
of sampled sequences, and fitness of the sampled sequences. If latent
cells are sampled, then an alignment of the sampled latent cells will
also be returned.

## Details

Also note that some of the inputs are expected to be in units of
*generations* and some are expected to be in units of *days*. These
choices were made based on what empirical estimates are most often
estimated present in the literature. We have highlighted in the
parameter descriptions which inputs are which.

## Examples

``` r
if (FALSE) { # \dontrun{
run_wavess(
  define_growth_curve(n_gen = 50),
  define_sampling_scheme(
    sampling_frequency_active = 10,
    sampling_frequency_latent = 10, n_days = 50
  ),
  rep("ATCG", 10)
)
} # }
```
