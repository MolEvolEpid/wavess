## code to prepare `hxb2_founder` dataset goes here

hxb2_founder <- ape::read.FASTA('data-raw/gp120_hxb2_founder_aligned.fasta')

usethis::use_data(hxb2_founder, overwrite = TRUE)
