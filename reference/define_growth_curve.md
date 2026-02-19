# Define active infected cell growth curve

Get per-generation logistic growth curve for infected active cells. Note
that the simulation will only be run until the final sampling time
within n_gens. See
[`define_sampling_scheme()`](https://molevolepid.github.io/wavess/reference/define_sampling_scheme.md).

## Usage

``` r
define_growth_curve(
  n_gens = 5000,
  n0 = 10,
  carry_cap = 2000,
  max_growth_rate = 0.3
)
```

## Arguments

- n_gens:

  number of generations (default: 5000)

- n0:

  starting infected cell population size (default: 10)

- carry_cap:

  carrying capacity for number of infected cells (default: 2000)

- max_growth_rate:

  maximum infected cell population growth rate (default: 0.3)

## Value

tibble with two columns: day and active cell count

## Examples

``` r
define_growth_curve()
#> # A tibble: 5,001 × 2
#>    generation active_cell_count
#>         <dbl>             <dbl>
#>  1          0                10
#>  2          1                13
#>  3          2                17
#>  4          3                22
#>  5          4                29
#>  6          5                37
#>  7          6                48
#>  8          7                62
#>  9          8                80
#> 10          9               103
#> # ℹ 4,991 more rows
```
