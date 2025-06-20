
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wavess

**Within-host agent-based virus evolution sequence simulator.**

<!-- <a href='https://github.com/MolEvolEpid/wavess/'><img src='man/figures/logo.png' align="right" height="120" /></a> -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/MolEvolEpid/wavess/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MolEvolEpid/wavess/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/MolEvolEpid/wavess/branch/main/graph/badge.svg)](https://app.codecov.io/gh/MolEvolEpid/wavess?branch=main)
<!-- badges: end -->

The goal of wavess is to simulate within-host virus sequence evolution
optionally including recombination, a latent infected cell reservoir,
and three types of selection (conserved sites, comparison to a fit
sequence, and antibody-mediated immunity). The package also provides
functions to pre-process data for input into the simulator, as well as
post-processing functions to analyze the simulation output. The
post-processing functions can also be used on real data. The default
settings for the simulator assume that the sequences being simulated are
HIV gp120.

Website: [molevolepid.github.io/wavess/](https://molevolepid.github.io/wavess/)

## Installation

You can install the development version of wavess from
[GitHub](https://github.com/MolEvolEpid/wavess) with:

``` r
# install.packages("devtools")
devtools::install_github("MolEvolEpid/wavess")
```

## Vignettes

Please see the following vignettes for examples:

- Preparing input data: `vignette("prepare_input_data")`
- Running wavess: `vignette("run_wavess")`
- Analyzing the output: `vignette("analyze_output")`
- Running the python script: `vignette("python")`
