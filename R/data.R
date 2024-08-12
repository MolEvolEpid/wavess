#' HIV1 ENV alignment
#'
#' The 2021 HIV1 ENV FLT web alignment from the LANL HIV database
#'
#' @format ## `hiv_env_flt_2021`
#' An `ape` DNAbin object with 6,741 DNA sequences of length 3993
#' @source <https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html>
"hiv_env_flt_2021"

#' gp120 HXB2 and founder alignment
#'
#' Alignment of gp120 HXB2 and DEMB11US006 (founder)
#'
#' @format ## `hxb2_founder`
#' An `ape` DNAbin object with 2 DNA sequences of length 1563 including the
#' env gp120 sequence for hxb2 and DEMB11US006 (founder)
#' @source
#' HXB2 (1st sequence):
#' <https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html>
#' Founder:
#' <https://www.sciencedirect.com/science/article/pii/S0022175914000143?via%3Dihub>
#' <https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=KC473833>
"hxb2_founder"

#' env features
#'
#' Binding, contact, and neutralization features in HIV env gene
#'
#' @format ## `env_features`
#' A tibble of HIV env binding, contact, and neutralization sites curated by the
#' LANL HIV database. The full Env feature dataset was subsetted to include only
#' binding, contact, and neutralization sites. The most relevant column for us
#' is the "position" column. Other columns are included in case it is useful for
#' specific use cases.
#' @source
#' https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/download_db.comp
"env_features"
