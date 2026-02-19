# Calculate sequence-based diversity and divergence

See
[`vignette('analyze_output')`](https://molevolepid.github.io/wavess/articles/analyze_output.md)
for more details.

## Usage

``` r
calc_div_metrics(aln, founder, gen)
```

## Arguments

- aln:

  DNA alignment in
  [ape::DNAbin](https://rdrr.io/pkg/ape/man/DNAbin.html) format

- founder:

  Name of the founder sequence in the alignment

- gen:

  Vector that indicates the generation of each sequence in the
  alignment, assumed to be in the same order as the alignment

## Value

tibble including mean sequence-based diversity and divergence for each
generation

## Examples

``` r
# This example is somewhat contrived, but it shows how it works.
calc_div_metrics(
  hxb2_cons_founder,
  "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
  c(1, 2, 2)
)
#> # A tibble: 2 × 3
#>     gen diversity divergence
#>   <dbl>     <dbl>      <dbl>
#> 1     1  NaN          0     
#> 2     2    0.0660     0.0542
```
