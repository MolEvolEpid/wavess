# Calculate nucleotide substitution probabilities

For this, you should probably use a set of within-host sequences for the
genome region of interest.

## Usage

``` r
estimate_q(aln, tr = NULL, model = "GTR+R(4)+I", rearrangement = "none")
```

## Arguments

- aln:

  An alignment of class ape::DNAbin

- tr:

  Starting tree

- model:

  A string providing model (e.g. "GTR+G(4)+I")

- rearrangement:

  Type of tree rearrangements to perform, one of "none", "NNI",
  "stochastic" or "ratchet" (default: "none")

## Value

Matrix of nucleotide substitution probabilities. Columns are from and
rows are to.

## Examples

``` r
# NOTE: This is just an example.
estimate_q(hiv_env_flt_2022)
#> optimize edge weights:  -13347 --> -13272.97 
#> optimize rate matrix:  -13272.97 --> -12963.97 
#> optimize invariant sites:  -12963.97 --> -12460.48 
#> optimize free rate parameters:  -12460.48 --> -10603 
#> optimize edge weights:  -10603 --> -10590.63 
#> optimize rate matrix:  -10590.63 --> -10581.46 
#> optimize invariant sites:  -10581.46 --> -10581.46 
#> optimize free rate parameters:  -10581.46 --> -10577.81 
#> optimize edge weights:  -10577.81 --> -10577.48 
#> optimize rate matrix:  -10577.48 --> -10577.26 
#> optimize invariant sites:  -10577.26 --> -10577.26 
#> optimize free rate parameters:  -10577.26 --> -10575.89 
#> optimize edge weights:  -10575.89 --> -10575.83 
#> optimize rate matrix:  -10575.83 --> -10575.79 
#> optimize invariant sites:  -10575.79 --> -10575.79 
#> optimize free rate parameters:  -10575.79 --> -10575.18 
#> optimize edge weights:  -10575.18 --> -10575.16 
#> optimize rate matrix:  -10575.16 --> -10575.14 
#> optimize invariant sites:  -10575.14 --> -10575.14 
#> optimize free rate parameters:  -10575.14 --> -10574.92 
#> optimize edge weights:  -10574.92 --> -10574.91 
#> optimize rate matrix:  -10574.91 --> -10574.91 
#> optimize invariant sites:  -10574.91 --> -10574.91 
#> optimize free rate parameters:  -10574.91 --> -10574.84 
#> optimize edge weights:  -10574.84 --> -10574.84 
#> optimize rate matrix:  -10574.84 --> -10574.83 
#> optimize invariant sites:  -10574.83 --> -10574.83 
#> optimize free rate parameters:  -10574.83 --> -10574.82 
#> optimize edge weights:  -10574.82 --> -10574.81 
#> optimize rate matrix:  -10574.81 --> -10574.81 
#> optimize invariant sites:  -10574.81 --> -10574.81 
#> optimize free rate parameters:  -10574.81 --> -10574.81 
#> optimize edge weights:  -10574.81 --> -10574.81 
#> optimize rate matrix:  -10574.81 --> -10574.81 
#> optimize invariant sites:  -10574.81 --> -10574.81 
#> optimize free rate parameters:  -10574.81 --> -10574.8 
#> optimize edge weights:  -10574.8 --> -10574.8 
#> optimize rate matrix:  -10574.8 --> -10574.8 
#> optimize invariant sites:  -10574.8 --> -10574.8 
#> optimize free rate parameters:  -10574.8 --> -10574.8 
#> optimize edge weights:  -10574.8 --> -10574.8 
#> optimize rate matrix:  -10574.8 --> -10574.8 
#> optimize invariant sites:  -10574.8 --> -10574.8 
#> optimize free rate parameters:  -10574.8 --> -10574.8 
#> optimize edge weights:  -10574.8 --> -10574.8 
#>            A          C          G          T
#> A -1.5940889  0.3894774  1.0395679  0.1650437
#> C  0.7692561 -2.0446348  0.1948555  1.0805232
#> G  1.5490726  0.1470087 -1.9302495  0.2341682
#> T  0.2459336  0.8152003  0.2341682 -1.2953020
```
