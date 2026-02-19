# Define sampling scheme

Define which days to sample sequences, and how many sequences to sample

## Usage

``` r
define_sampling_scheme(
  sampling_frequency_active = 365,
  max_samp_active = 20,
  sampling_frequency_latent = 365,
  max_samp_latent = 20,
  n_days = 3650
)
```

## Arguments

- sampling_frequency_active:

  frequency in days at which to record sequences from active cells (and
  counts) (default: 365 days)

- max_samp_active:

  maximum number of cells (and thus sequences) to sample from active
  cells in a given day (default: 20 sequences)

- sampling_frequency_latent:

  frequency in days at which to record sequences from latent cells
  (default: 365 days)

- max_samp_latent:

  maximum number of cells (and thus sequences) to sample from latent
  cells in a given day (default: 20 sequences)

- n_days:

  day to end sampling (default: 3650)

## Value

input growth curve tibble with one additional column (`n_sample_active`)
containing the number of sequences from active cells to samples

## Examples

``` r
define_sampling_scheme()
#> # A tibble: 11 × 3
#>      day n_sample_active n_sample_latent
#>    <int>           <dbl>           <dbl>
#>  1     0              20              20
#>  2   365              20              20
#>  3   730              20              20
#>  4  1095              20              20
#>  5  1460              20              20
#>  6  1825              20              20
#>  7  2190              20              20
#>  8  2555              20              20
#>  9  2920              20              20
#> 10  3285              20              20
#> 11  3650              20              20
```
