# Calculate Q matrix from nucleotide substitution rates

Method:

- If needed, convert rates per day to rates per generation

- Convert rates per generation to probabilities per generation

- Make diagonal such that rows of probability matrix p sum to 1

- Convert to q matrix, assuming a mutation rate of mut_rate
  mutations/site/generation, by solving for q in the equation p =
  exp(q\*mut_rate)

## Usage

``` r
calc_q_from_rates(rates, mut_rate, generation_time = NULL)
```

## Arguments

- rates:

  4x4 matrix of individual nucleotide substitution rates for each
  potential substitution (per-site per-generation or per-site per-day)
  with the following row and column names: A,C,G,T.

- mut_rate:

  Overall per-site per-generation mutation rate.

- generation_time:

  Generation time (if nucleotide substitution rates are in days rather
  than generations; default: NULL, assumes that rate matrix is
  per-generation)

## Value

Nucletoide substitution rate matrix Q

## Examples

``` r
calc_q_from_rates(hiv_mut_rates, 2.4e-5, 1.2)
#>            A          C          G          T
#> A -1.3935602  0.1650215  1.1002145  0.1283242
#> C  0.9167813 -3.2089352  0.0915972  2.2005567
#> G  2.9338881  0.0182889 -3.3189236  0.3667466
#> T  0.5499598  1.8338005  0.5501132 -2.9338736
```
