#' HIV1 ENV alignment
#'
#' The 2021 HIV1 ENV FLT web alignment from the LANL HIV database
#'
#' @format ## `hiv_env_flt_2021`
#' An `ape` DNAbin object with 10 DNA sequences of length 3993.
#' If you would like to download the entire set of 6,741 sequences, please
#' visit the link below.
#' @source <https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html>
"hiv_env_flt_2021"

#' HXB2, consensus, and founder alignment
#'
#' Full-genome alignment of HXB2, consensus sequence, and DEMB11US006 (founder)
#'
#' @format ## `hxb2_cons_founder` An `ape` DNAbin object with 3 DNA sequences of
#'   length 1563 including the full-genome sequence for hxb2, the consensus
#'   sequence, and DEMB11US006 (founder)
#' @source HXB2 and consensus (1st 2 sequences):
#' <https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html>
#' Founder:
#' <https://www.sciencedirect.com/science/article/pii/S0022175914000143?via%3Dihub>
#' <https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=KC473833>
"hxb2_cons_founder"

#' Example founder conserved sites
#'
#' Conserved sites (indexed at 1) for the gp120 DEMB11US006 sequence, which is
#' used as an example founder sequence throughout the package. These sites were
#' identified using the `identify_conserved_sites()` function with the full
#' HIV ENV filtered alignment from the LANL HIV database.
#' See `hiv_env_flt_2021` for more details and the link to the entire alignment.
#'
#' @format ## `conserved_sites` Vector of conserved nucleotide sites for the
#' gp120 gene of DEMB11US006 (founder)
#' @source Full alignment used:
#' <https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html>
#' Founder sequence:
#' <https://www.sciencedirect.com/science/article/pii/S0022175914000143?via%3Dihub>
#' <https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=KC473833>
"conserved_sites"

#' ENV features
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
