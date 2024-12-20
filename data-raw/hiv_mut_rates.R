## generate hiv q matrix

# per-day mutation rates based on approximately neutral sites from Fig 1C in manuscript:
# Fabio Zanini, Vadim Puller, Johanna Brodin, Jan Albert, Richard A. Neher,
# In vivo mutation rates and the landscape of fitness costs of HIV-1, Virus Evolution,
# Volume 3, Issue 1, January 2017, vex003, https://doi.org/10.1093/ve/vex003
hiv_mut_rates <- as.matrix(data.frame(
  A = c(0, 5e-6, 1.6e-5, 3e-6),
  C = c(9e-7, 0, 1e-7, 1e-5),
  G = c(6e-6, 5e-7, 0, 3e-6),
  T = c(7e-7, 1.2e-5, 2e-6, 0)
))
usethis::use_data(hiv_mut_rates, overwrite = TRUE)


