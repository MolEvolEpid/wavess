
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wavess

**Within-host agent-based viral evolution sequence simulator.**

<!-- badges: start -->

[![R-CMD-check](https://github.com/MolEvolEpid/wavess/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MolEvolEpid/wavess/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of wavess is to simulate within-host viral sequence evolution
optionally including recombination, a latent infected cell reservoir,
and three types of selection (conserved sites, comparison to a fit
sequence, and antibody-mediated immunity). The package also provides
functions to pre-process data for input into the simulator, as well as
post-processing functions to analyze the simulation output. Theh
post-processing functions can also be used on real data. The default
settings for the simulator assume that the sequences being simulated are
HIV gp120.

## Installation

You can install the development version of wavess from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MolEvolEpid/wavess")
```

## Vignettes

Please see the vignettes following vignettes for examples:

- Preparing input data: `vignette("prepare_input_data")`
- Running wavess: `vignette("run_wavess")`
- Analyzing the output: `vignett("analyze_output")`
