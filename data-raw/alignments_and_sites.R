## code to prepare `hiv_env_flt_2022` dataset
# library(wavess)

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
  extract_seqs(hcf_gp120,"B.FR.83.HXB2_LAI_IIIB_BRU.K03455")$founder
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
                           founder_aln = hcf_gp120) |>
  dplyr::filter(conserved == "Yes") |>
  dplyr::pull(founder_pos)

usethis::use_data(conserved_sites, compress = "xz", overwrite = TRUE)
