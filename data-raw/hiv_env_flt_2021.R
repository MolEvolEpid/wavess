## code to prepare `hiv_env_flt_2021` dataset goes here

hiv_env_flt_2021 <- as.matrix(ape::read.FASTA("data-raw/HIV1_FLT_2021_env_DNA.fasta"))[1:10, ]

usethis::use_data(hiv_env_flt_2021, compress = "xz", overwrite = TRUE)
