## code to prepare `hxb2_cons_founder` dataset goes here

hxb2_cons_founder <- as.matrix(ape::read.FASTA("data-raw/hxb2_cons_founder_aligned.fasta"))

usethis::use_data(hxb2_cons_founder, compress = "xz", overwrite = TRUE)
