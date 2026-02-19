# Get epitope frequencies

Get epitope frequencies

## Usage

``` r
get_epitope_frequencies(epitope_positions)
```

## Arguments

- epitope_positions:

  numeric vector of **amino acid** positions of identified epitopes
  (e.g. the env features from the [LANL HIV
  database](https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/download_db.comp))

## Value

tibble with the following columns:

- `aa_position`: amino acid position

- `n_features`: number of features at that position

- `epitope_probability`: estimated probability of an epitope at that
  position

## Examples

``` r
get_epitope_frequencies(env_features$Position - 1) # subtract 1 to index at 0
#> # A tibble: 480 × 3
#>    aa_position n_features epitope_probability
#>          <dbl>      <dbl>               <dbl>
#>  1          30          1            0.000331
#>  2          31          0            0       
#>  3          32          0            0       
#>  4          33          0            0       
#>  5          34          0            0       
#>  6          35          0            0       
#>  7          36          0            0       
#>  8          37          0            0       
#>  9          38          0            0       
#> 10          39          0            0       
#> # ℹ 470 more rows
```
