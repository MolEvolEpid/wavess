# Sample epitopes

Sample epitopes based on epitope probabilities. Note that the positions
returned assume that the start of the amino acid sequence is also the
start of the founder sequence in the simulation. We also assume that
there are no frameshift mutations in the founder sequence.

## Usage

``` r
sample_epitopes(
  epitope_probabilities,
  start_aa_pos = 0,
  end_aa_pos = NULL,
  num_epitopes = 10,
  aa_epitope_length = 10,
  max_fit_cost = 0.3,
  max_resamples = 100,
  ref_founder_map = NULL
)
```

## Arguments

- epitope_probabilities:

  Epitope probability tibble as output by
  [`get_epitope_frequencies()`](https://molevolepid.github.io/wavess/reference/get_epitope_frequencies.md),
  including columns `aa_position` and `epitope_probability`.
  `aa_position` should be indexed at 0

- start_aa_pos:

  Starting amino acid position to consider for epitopes, indexed at 0
  (default: 0, i.e. the first position)

- end_aa_pos:

  Ending amino acid position to consider for epitopes, indexed at 0
  (default: NULL, i.e. through the final position in
  `epitope_probabilities$aa_position`)

- num_epitopes:

  Number of epitopes to sample

- aa_epitope_length:

  Amino acid epitope length

- max_fit_cost:

  Maximum fitness cost of an epitope, must be in the range \[0,1) where
  0 indicates no cost. 1, which indicates no ability to survive, is not
  allowed (default: 0.3) **note that the model output is very sensitive
  to this parameter**

- max_resamples:

  Maximum number of resampling events to attempt; this is to prevent an
  infinite loop (default: 100)

- ref_founder_map:

  Output from
  [`map_ref_founder()`](https://molevolepid.github.io/wavess/reference/map_ref_founder.md),
  including *nucleotide* reference and founder positions (`ref_pos` and
  `founder_pos`). **NOTE:** The reference positions here, if they were
  converted to amino acid positions, are expected to match with the
  reference positions in `epitope_probabilities`. Further, we assume
  that the founder indices align with the founder sequence positions to
  be used in the simulation (default: NULL)

## Value

tibble with the `num_epitopes` rows and the following columns:

- `epi_start_nt`: nucleotide epitope start position

- `epi_end_nt`: nucleotide epitope end position

- `max_fitness_cost`: maximum fitness cost for that epitope

## Examples

``` r
sample_epitopes(get_epitope_frequencies(env_features$Position - 1))
#> 2 resamples required
#> # A tibble: 10 × 3
#>    epi_start_nt epi_end_nt max_fitness_cost
#>           <dbl>      <dbl>            <dbl>
#>  1          273        303             0.03
#>  2          486        516             0.06
#>  3          909        939             0.09
#>  4          384        414             0.12
#>  5          456        486             0.15
#>  6          702        732             0.18
#>  7         1173       1203             0.21
#>  8         1251       1281             0.24
#>  9         1092       1122             0.27
#> 10          984       1014             0.3 
```
