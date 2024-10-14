#' Create python virtual environment and dependencies in conda or virtualenv
#' environment
#'
#' Modified from the [spacyr::spacy_install()] function:
#' https://github.com/quanteda/spacyr/blob/master/R/spacy_install.R
#'
#' The dependencies are python, some base python packages (copy, random,
#' collections), and numpy.
#'
#' @param ask logical; ask whether to proceed during the installation. By
#'   default, questions are only asked in interactive sessions.
#' @param force ignore if the dependencies are already present and install it
#'   anyway.
#'
#' @details The function checks whether a suitable installation of Python is
#'   present on the system and installs one via [reticulate::install_python()]
#'   otherwise. It then creates a virtual environment with the necessary
#'   packages in the default location chosen by [reticulate::virtualenv_root()].
#'
#'   If you want to install a different version of Python than the default, you
#'   should call [reticulate::install_python()] directly. If you want to create
#'   or use a different virtual environment, you can use, e.g.,
#'   `Sys.setenv(WAVESS_PYTHON = "path/to/directory")`. If you want to install a
#'   different version of a particular package, then use, e.g.
#'   reticulate::virtualenv_install(packages = c("numpy==VERSION_NUM")) instead.
#'
#'
#' @examples
#' \dontrun{
#' # install dependencies
#' create_python_venv()
#'
#' # update dependencies
#' create_python_venv(force = TRUE)
#' }
#'
#' @export
create_python_venv <- function(ask = interactive(),
                               force = FALSE) {
  if (!reticulate::virtualenv_exists(Sys.getenv("WAVESS_PYTHON",
    unset = "r-wavess"
  ))) {
    # this has turned out to be the easiest way to test if a suitable Python
    # version is present. All other methods load Python, which creates
    # some headache.
    t <- try(
      reticulate::virtualenv_create(Sys.getenv("WAVESS_PYTHON",
        unset = "r-wavess"
      )),
      silent = TRUE
    )
    if (methods::is(t, "try-error")) {
      permission <- TRUE
      if (ask) {
        permission <- utils::askYesNo(paste0(
          "No suitable Python installation was found on your system. ",
          "Do you want to run `reticulate::install_python()` to install it?"
        ))
      }

      if (permission) {
        if (utils::packageVersion("reticulate") < "1.19") {
          stop(
            "Your version or reticulate is too old for this action. ",
            "Please update"
          )
        }
        python <- reticulate::install_python()
        reticulate::virtualenv_create(
          Sys.getenv("WAVESS_PYTHON",
            unset = "r-wavess"
          ),
          python = python
        )
      } else {
        stop("Aborted by user")
      }
    }
  }
  reticulate::use_virtualenv(Sys.getenv("WAVESS_PYTHON", unset = "r-wavess"))

  dependencies <- c("numpy")

  installed <- c()
  for (dep in dependencies) {
    if (py_check_installed(dep) && !force) {
      installed <- c(installed, dep)
    } else {
      reticulate::py_install(dep, envname = Sys.getenv("WAVESS_PYTHON",
        unset = "r-wavess"
      ))
      message(
        "Installation of ", dep, " version ",
        py_check_version(dep, envname = Sys.getenv("WAVESS_PYTHON",
          unset = "r-wavess"
        )),
        " complete."
      )
    }
  }

  if (length(installed) > 0) {
    warning(
      "Skipped installation of the following packages: ",
      installed,
      " Use `force` to force installation or update."
    )
  }
  invisible(NULL)
}

#' Uninstall the wavess environment
#'
#' Removes the virtual environment created by create_python_venv()
#'
#' @param confirm logical; confirm before uninstalling wavess virtual
#'   environment?
#'
#' @export
remove_python_venv <- function(confirm = interactive()) {
  reticulate::virtualenv_remove(Sys.getenv("WAVESS_PYTHON", unset = "r-wavess"),
    confirm = confirm
  )

  message("Deinstallation complete.")
  invisible(NULL)
}

#' Check if python dependency is installed
#'
#' @param x dependency to check
#'
#' @return boolean whether dependency is installed
#' @noRd
py_check_installed <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  return(x %in% trimws(reticulate::py_list_packages(
    Sys.getenv("WAVESS_PYTHON", unset = "r-wavess")
  )$package))
}

#' Check python package version
#'
#' @param package Package to check version of
#' @param ... Extra stuff for reticulate::py_list_packages
#'
#' @return package version
#' @noRd
py_check_version <- function(package, ...) {
  packages <- reticulate::py_list_packages(...)
  packages$version[packages$package == package]
}


#' Use wavess python virtual environment
#'
#' @return functions in agents
#' @noRd
use_python_venv <- function() {
  if (!reticulate::virtualenv_exists(
    Sys.getenv("WAVESS_PYTHON", unset = "r-wavess")
  )) {
    stop(
      "No wavess environment found. ",
      "Use `install_python_dependencies()` to get started."
    )
  }

  if (!"numpy" %in% reticulate::py_list_packages(
    Sys.getenv("WAVESS_PYTHON", unset = "r-wavess")
  )$package) {
    stop(
      "numpy was not found in your environment. ",
      "Use `install_python_dependencies()`",
      "to get started."
    )
  }

  reticulate::use_virtualenv(Sys.getenv("WAVESS_PYTHON", unset = "r-wavess"))
  return(reticulate::import_from_path("agents",
    system.file("python", package = "wavess"),
    convert = FALSE
  ))
}
