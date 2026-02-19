# Create python virtual environment and dependencies in conda or virtualenv environment

Modified from the `spacyr::spacy_install()` function:
https://github.com/quanteda/spacyr/blob/master/R/spacy_install.R

## Usage

``` r
create_python_venv(ask = interactive(), force = FALSE)
```

## Arguments

- ask:

  logical; ask whether to proceed during the installation. By default,
  questions are only asked in interactive sessions.

- force:

  ignore if the dependencies are already present and install it anyway.

## Details

The dependencies are python, some base python packages (copy, random,
collections), numpy, scipy, and pandas.

The function checks whether a suitable installation of Python is present
on the system and installs one via
[`reticulate::install_python()`](https://rstudio.github.io/reticulate/reference/install_python.html)
otherwise. It then creates a virtual environment with the necessary
packages in the default location chosen by
[`reticulate::virtualenv_root()`](https://rstudio.github.io/reticulate/reference/virtualenv-tools.html).

If you want to install a different version of Python than the default,
you should call
[`reticulate::install_python()`](https://rstudio.github.io/reticulate/reference/install_python.html)
directly. If you want to create or use a different virtual environment,
you can use, e.g., `Sys.setenv(WAVESS_PYTHON = "path/to/directory")`. If
you want to install a different version of a particular package, then
use, e.g. reticulate::virtualenv_install(packages =
c("numpy==VERSION_NUM")) instead.

## Examples

``` r
if (FALSE) { # \dontrun{
# install dependencies
create_python_venv()

# update dependencies
create_python_venv(force = TRUE)
} # }
```
