# wavess

**Within-host adaptive virus evolution sequence simulator.**

The goal of wavess is to simulate within-host virus sequence evolution
optionally including recombination, a latent infected cell reservoir,
and three types of selection (conserved sites, comparison to a fit
sequence, antibody-mediated immunity, and cytotoxic T-cell mediated
immunity). The package also provides functions to pre-process data for
input into the simulator, as well as post-processing functions to
analyze the simulation output. The post-processing functions can also be
used on real data. The default settings for the simulator assume that
the sequences being simulated are HIV gp120.

Website:
[molevolepid.github.io/wavess/](https://molevolepid.github.io/wavess/)

Manuscript: [wavess: An R package for simulation of adaptive within-host
virus sequence
evolution](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013437)

\*\*Note that another manuscript is in prep describing the CTL response
and variable recombination rate. Release v1.1 includes the variable
recombination rate. Upcoming release v1.2 will include the CTL response.
Installing as described below includes the CTL response implementation.

## Installation

You can install the development version of wavess from
[GitHub](https://github.com/MolEvolEpid/wavess) with:

``` r
# install.packages("devtools")
devtools::install_github("MolEvolEpid/wavess")
```

This will not install the vignettes locally, but you can still see the
vignettes on the [website](https://molevolepid.github.io/wavess/). If
you’d like to install the vignettes locally, you can use the command:

``` r
# install.packages("devtools")
devtools::install_github("MolEvolEpid/wavess", build_vignettes = TRUE, force = TRUE)
```

## Vignettes

Please go to the [website](https://molevolepid.github.io/wavess/) and
see the following vignettes for examples:

- Preparing input data:
  [`vignette("prepare_input_data")`](https://molevolepid.github.io/wavess/articles/prepare_input_data.md)
- Running wavess:
  [`vignette("run_wavess")`](https://molevolepid.github.io/wavess/articles/run_wavess.md)
- Analyzing the output:
  [`vignette("analyze_output")`](https://molevolepid.github.io/wavess/articles/analyze_output.md)
- Running the python script:
  [`vignette("python")`](https://molevolepid.github.io/wavess/articles/python.md)

## Citation

``` r
citation("wavess")
```

## Contributing

You are of course welcome to fork the package for your own use. If you
would like for your changes to be incorporated into the latest version
of `wavess`, please open an
[issue](https://github.com/MolEvolEpid/wavess/issues) related to your
proposed change, create a [pull
request](https://github.com/MolEvolEpid/wavess/pulls) with the
implemented change, or reach out to <tkl@lanl.gov>,
<nsambaturu@binghamton.edu>, or <zenalapp@lanl.gov>.
