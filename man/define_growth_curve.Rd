% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_pop_samp.R
\name{define_growth_curve}
\alias{define_growth_curve}
\title{Define active infected cell growth curve}
\usage{
define_growth_curve(
  n_gens = 5000,
  n0 = 10,
  carry_cap = 2000,
  max_growth_rate = 0.3
)
}
\arguments{
\item{n_gens}{number of generations (default: 5000)}

\item{n0}{starting infected cell population size (default: 10)}

\item{carry_cap}{carrying capacity for number of infected cells (default: 2000)}

\item{max_growth_rate}{maximum infected cell population growth rate
(default: 0.3)}
}
\value{
tibble with two columns: day and active cell count
}
\description{
Get per-generation logistic growth curve for infected active cells.
Note that the simulation will only be run until the final sampling time
within n_gens. See \code{\link[=define_sampling_scheme]{define_sampling_scheme()}}.
}
\examples{
define_growth_curve()
}
