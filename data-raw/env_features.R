## code to prepare `env_features` dataset goes here

env_features <- readr::read_tsv("data-raw/env_features_2024-10-01.tsv") |>
  # filter to include only binding, contacts, and neutralization
  # don't include Env feature, immunotherapy, other feature, signature
  dplyr::filter(`Feature type` %in% c("binding", "contacts", "neutralization") &
    # only include features in gp120
    grepl("gp120", `Env feature(s)`))

usethis::use_data(env_features, compress = "xz", overwrite = TRUE)
