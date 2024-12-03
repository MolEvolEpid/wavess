## generate hiv q matrix

# rates based on approximately neutral sites from Fig 1C in manuscript:
# Fabio Zanini, Vadim Puller, Johanna Brodin, Jan Albert, Richard A. Neher,
# In vivo mutation rates and the landscape of fitness costs of HIV-1, Virus Evolution,
# Volume 3, Issue 1, January 2017, vex003, https://doi.org/10.1093/ve/vex003
# Method:
# - convert rates per day to rates per generation
# - convert rates per generation to probabilities per generation
# - make diagonal such that rows of probability matrix p sum to 1
# - convert to q matrix, assuming a mutation rate t of 3.5e-5 mutations/site/generation,
#   by solving for q in the equation p = exp(qt)
rates_per_day <- as.matrix(data.frame(
  A = c(0, 5e-6, 1.6e-5, 3e-6),
  C = c(9e-7, 0, 1e-7, 1e-5),
  G = c(6e-6, 5e-7, 0, 3e-6),
  T = c(7e-7, 1.2e-5, 2e-6, 0)
))
rates_per_gen <- rates_per_day * 1.2
probs_per_gen <- rate_to_probability(rates_per_gen)
diag(probs_per_gen) <- 1 - rowSums(probs_per_gen)
hiv_q_mat <- expm::logm(probs_per_gen) / 3.5e-5
rownames(hiv_q_mat) <- colnames(hiv_q_mat) <- c("A", "C", "G", "T")
usethis::use_data(hiv_q_mat, overwrite = TRUE)
