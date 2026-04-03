# Calculate tree summary statistics

See
[`vignette('analyze_output')`](https://molevolepid.github.io/wavess/articles/analyze_output.md)
for more details.

## Usage

``` r
calc_tr_stats(tr, timepoints, bl_thresh = 1e-08)
```

## Arguments

- tr:

  Rooted phylogeny

- timepoints:

  Factor vector of time points named by tree tip label. The levels
  should be ordered correctly (most often by sampling time). All tip
  lables must be included in this vector. If you want to exclude certain
  tips, you must drop them from the tree prior to using this function.

- bl_thresh:

  Branch length threshold under which branches are collapsed to prior to
  calculations (default: 1e-08). This is the `tol` argument in
  [`ape::di2multi()`](https://rdrr.io/pkg/ape/man/multi2di.html).

## Value

Tibble including tree summary statistics:

- Mean leaf depth (normalized Sackin index)

- Mean branch length

- Mean internal branch length

- Mean external branch length

- Mean divergence (mean per-generation root-to-tip distance)

- Mean diversity (mean per-generation tip-to-tip distance)

- Divergence slope across timepoints

- Diversity slope across timepoints

- Mean number of lineages through time

## Examples

``` r
tr <- ape::rtree(100)
times <- factor(sample(3, 100, replace = TRUE), levels = 1:3)
names(times) <- tr$tip.label
calc_tr_stats(tr, times)
#> # A tibble: 9 × 2
#>   stat_name        stat_value
#>   <chr>                 <dbl>
#> 1 mean_leaf_depth     8.43   
#> 2 mean_bl             0.527  
#> 3 mean_int_bl         0.509  
#> 4 mean_ext_bl         0.544  
#> 5 mean_divergence     4.22   
#> 6 mean_diversity      7.00   
#> 7 divergence_slope   -0.00326
#> 8 diversity_slope    -0.0413 
#> 9 mean_ltt           22.5    
```
