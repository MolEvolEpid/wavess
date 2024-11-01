## code to prepare `hiv_env_flt_2022` dataset
library(wavess)

# alignment with hxb2, consensus, and founder
hxb2_cons_founder <- as.matrix(
  ape::read.FASTA("data-raw/hxb2_cons_founder_aligned.fasta")
)

usethis::use_data(hxb2_cons_founder, compress = "xz", overwrite = TRUE)

# filtered alignment from LANL HIV database
hiv_env_flt_2022_complete <- as.matrix(
  ape::read.FASTA("data-raw/HIV1_FLT_2022_env_DNA.fasta")
)

# subset for examples
hiv_env_flt_2022 <- hiv_env_flt_2022_complete[1:10, ]
usethis::use_data(hiv_env_flt_2022, compress = "xz", overwrite = TRUE)

# slice out gp120 for each alignment

hcf_gp120_start <- which(cumsum(as.character(
  hxb2_cons_founder["B.FR.83.HXB2_LAI_IIIB_BRU.K03455", ]
) != "-") == 6225)[1]
hcf_gp120_end <- which(cumsum(as.character(
  hxb2_cons_founder["B.FR.83.HXB2_LAI_IIIB_BRU.K03455", ]
) != "-") == 7757)[1]
hcf_gp120 <- slice_aln(hxb2_cons_founder, hcf_gp120_start, hcf_gp120_end)

len_hxb2_gp120 <- nchar(
  extract_seqs(hcf_gp120, "B.FR.83.HXB2_LAI_IIIB_BRU.K03455")$founder
)
flt_gp120_start <- 1
flt_gp120_end <- which(cumsum(as.character(hiv_env_flt_2022_complete[1, ]) !=
  "-") == len_hxb2_gp120)[1]
# subset to only gp120 section
flt_gp120 <- slice_aln(hiv_env_flt_2022_complete, 1, flt_gp120_end)

conserved_sites <-
  identify_conserved_sites(flt_gp120,
    founder = "B.US.2011.DEMB11US006.KC473833",
    ref = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
    founder_aln = hcf_gp120
  ) |>
  dplyr::filter(conserved == "Yes") |>
  dplyr::select(founder_pos, founder_base) |>
  tibble::deframe() |>
  toupper()

usethis::use_data(conserved_sites, compress = "xz", overwrite = TRUE)

## generate hiv q matrix

# rates from Fig 1C in manuscript:
# Fabio Zanini, Vadim Puller, Johanna Brodin, Jan Albert, Richard A. Neher,
# In vivo mutation rates and the landscape of fitness costs of HIV-1, Virus Evolution,
# Volume 3, Issue 1, January 2017, vex003, https://doi.org/10.1093/ve/vex003
# need to:
# - convert per-site per-day rates to per-site per-generation (1.2 days) rates
# - multiply by "to" nucleotide frequencies to get the Q matrix
# - divide by the overall substitution rate per site per generation to get the
#   unitless Q matrix
# TODO - THINK ABOUT WHETHER THIS IS THE RIGHT WAY TO DO THIS
nt_freqs <- ape::base.freq(hiv_env_flt_2022_complete)
names(nt_freqs) <- toupper(names(nt_freqs))
hiv_q_mat <- as.matrix(data.frame(
  A = c(0, 5e-6, 1.6e-5, 3e-6) * nt_freqs["A"],
  C = c(9e-7, 0, 1e-7, 1e-5) * nt_freqs["C"],
  G = c(6e-6, 5e-7, 0, 3e-6) * nt_freqs["G"],
  T = c(7e-7, 1.2e-5, 2e-6, 0) * nt_freqs["T"]
)) * 1.2 / 3.5e-5
diag(hiv_q_mat) <- -rowSums(hiv_q_mat)
rownames(hiv_q_mat) <- colnames(hiv_q_mat)

usethis::use_data(hiv_q_mat, overwrite = TRUE)
