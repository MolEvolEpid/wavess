#' Identify conserved sites
#'
#' Identify conserved sites relative to the founder sequence.
#' Note that the positions returned assume that the start of the alignment
#' is also the start of the founder sequence in the simulation.
#'
#' @param thresh Conserved site threshold. A position is considered to be
#' conserved if >thresh proportion of sequences in the alignment are the same
#' base (default: 0.99)
#'
#' @inheritParams find_consensus
#'
#' @return Tibble with the same columns as output in `find_consensus()`, plus
#' an additional column, `conserved`, that indicates whether or not the position
#' is conserved (Yes means conserved, No means not conserved, NA means the
#' conserved position is a gap ('-'))
#'
#' @export
#'
#' @examples
#' hiv_gp120_flt_2021 <- slice_aln(hiv_env_flt_2021, start = 1, end = 2517)
#' hxb2_cons_founder <- slice_aln(hxb2_cons_founder, start = 6225, end = 7757)
#' identify_conserved_sites(hiv_gp120_flt_2021, 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455')
#' identify_conserved_sites(hiv_gp120_flt_2021, 'B.US.2011.DEMB11US006.KC473833',
#' ref = 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', founder_aln = hxb2_cons_founder)
identify_conserved_sites <- function(aln, founder, thresh = 0.99,
                                     ref = NULL, founder_aln = NULL){
  check_identify_conserved_sites_inputs(aln, founder, thresh, ref, founder_aln)
  conserved_dat <- find_consensus(aln, founder, ref, founder_aln) |>
    # Remove reference positions that are gaps in >thresh sequences

    # dplyr::filter(!is.na(.data$consensus_base)) |>
    # Consider conserved sequences to be those that are identical in >thresh sequences
    # But if it's a gap, return the consensus base to be NA
    # because they would erroneously be called conserved in reference
    dplyr::mutate(conserved = dplyr::case_when(.data$consensus_base == '-' & .data$consensus_prop > thresh ~ NA,
                                               .data$consensus_prop > thresh ~ 'Yes',
                                               TRUE ~ 'No')) |>
    dplyr::select('founder_pos', 'founder_base', 'consensus_base', 'consensus_prop', 'conserved')

  return(conserved_dat)
}

#' Determine consensus sequence
#'
#' This function takes a very simple approach to determining the consensus
#' sequence. It simply finds the majority base. If there is a tie, it takes the
#' base that comes first alphabetically. For more nuanced consensus-making,
#' please check out the
#' [LANL HIV database consensus maker](https://www.hiv.lanl.gov/content/sequence/CONSENSUS/consensus.html)
#'
#' @param aln Alignment in `ape::DNAbin` format (e.g. read in with
#' `ape::read.dna()` or `ape::read.FASTA()`) that includes the founder sequence or
#' a reference sequence, plus a representative set of sequences for your genome
#' region of interest
#' @param founder  Name of the founder sequence in the input alignment or in the
#' optional founder-specific alignment (`founder_aln`)
#' @param ref Name of the reference sequence in the input alignment, required if
#' the founder sequence is not in `aln` (default: `NULL`)
#' @param founder_aln Alignment including the reference and founder sequences,
#' required if the founder is not present in `aln` (default: `NULL`)
#'
#' @return Tibble including the following columns:
#' - `founder_pos`: founder position
#' - `founder_base`: founder base
#' - `consensus_base`: consensus base
#' - `consensus_prop`: proportion of sequences that had that base at that position
#' When using a reference, `NA` in the consensus columns indicates that that
#' position was an insertion relative to the reference
find_consensus <- function(aln, founder, ref = NULL, founder_aln = NULL){
  check_find_consensus_inputs(aln, founder, ref, founder_aln)
  # ensure the alignment is in matrix form
  aln <- as.matrix(aln)

  if(is.null(ref)){
    ref <- founder
  }
  # get consensus base and proportion for each position
  consensus_info <- sapply(1:ncol(aln), function(x){
    tab <- table(as.character(aln[,x]))
    prop_tab <- prop.table(tab)
    c(consensus_base = names(tab)[which.max(prop_tab)],
      consensus_prop = max(prop_tab))
  })
  ref_consensus_dat <- tibble::tibble(ref = as.character(aln[ref,])[1,],
                                      consensus_base = consensus_info[1,],
                                      consensus_prop = consensus_info[2,]) |>
    get_seq_pos('ref') |>
    dplyr::rename(ref_base = 'ref') |>
    # Remove reference positions that are gaps
    # This might remove a small handful of bases that are conserved in other sequences,
    # but it would be complicated and low return to try to match these to other sequences
    dplyr::filter(!is.na(.data$ref_pos))
  if(is.null(founder_aln)){
    founder_consensus_dat <- ref_consensus_dat |>
      dplyr::rename(founder_base = 'ref_base', founder_pos = 'ref_pos')
  }else{
    # Consensus sites relative to founder
    ref_founder_map <- map_ref_founder(founder_aln, ref, founder)
    founder_consensus_dat <- ref_founder_map |>
      dplyr::select('founder_pos', 'founder_base') |>
      dplyr::filter(!is.na(.data$founder_pos)) |>
      dplyr::left_join(ref_founder_map |>
      dplyr::inner_join(ref_consensus_dat,
                        by = dplyr::join_by('ref_pos', 'ref_base')) |>
      # Remove gaps in founder
      dplyr::filter(!is.na(.data$founder_pos)),
      by = dplyr::join_by('founder_pos', 'founder_base')) |>
      dplyr::select(-c('alignment_pos', 'ref_pos', 'ref_base'))
  }
  founder_consensus_dat <- founder_consensus_dat |>
    dplyr::mutate(consensus_prop = as.numeric(.data$consensus_prop))
  return(founder_consensus_dat)
}


#' Map reference and founder sequence positions
#'
#' @param aln Alignment in `ape::DNAbin` format (e.g. read in with
#' `ape::read.dna()` or `ape::read.FASTA()`) that includes the reference sequence
#' and the founder sequence
#' @param ref Name of the reference sequence in the input alignment
#' @param founder Name of the founder sequence in the input alignment
#'
#' @return Tibble including columns for:
#' - Alignment position (`alignment_pos`)
#' - Reference base position (`ref_pos`)
#' - Founder base position (`founder_pos`)
#' - Reference base (`ref_base`)
#' - Founder base (`founder_base`)
#'
#' @export
#'
#' @examples
#' map_ref_founder(hxb2_cons_founder,
#' 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455',
#' 'B.US.2011.DEMB11US006.KC473833')
map_ref_founder <- function(aln, ref, founder){
  check_map_ref_founder_inputs(aln, ref, founder)
  # read in alignment and select reference and founder sequences
  aln <- as.list(aln)[c(ref, founder)]
  # rename sequences to make the next part easier
  names(aln) <- c('ref', 'founder')

  # do some transformations to get the alignment in a useful format
  aln |> as.matrix() |> as.character() |> t() |> tibble::as_tibble() |>
    # get alignment positions
    dplyr::mutate(alignment_pos = dplyr::row_number()) |>
    # get sequence positions for reference
    get_seq_pos('ref') |>
    # get sequence positions for founder
    get_seq_pos('founder') |>
    # remove leading and trailing positions relative to founder sequence
    dplyr::filter(cumsum(tidyr::replace_na(.data$founder_pos, 0)) > 0 &
             rev(cumsum(tidyr::replace_na(rev(.data$founder_pos), 0))) > 0) |>
    # select columns of interest
    dplyr::select('alignment_pos', 'ref_pos', 'founder_pos',
                  ref_base='ref', founder_base='founder')
}

#' Get sequence position
#'
#' Get sequence position from alignment data frame.
#' Currently considers everything a position that isn't '-'
#'
#' @param aln_df Alignment data frame including a column with the sequence.
#' @param col_name The column name of the sequence of interest.
#'
#' @return
#' Original data frame with an additional column (the column name prepended to
#' _pos) that indicates the position of that base in the sequence.n
get_seq_pos <- function(aln_df, col_name){
  check_get_seq_pos_inputs(aln_df, col_name)
  aln_df |>
    # determine where the gaps are
    dplyr::mutate(gap = .data[[col_name]] == '-') |>
    # get the sequence position
    dplyr::group_by(.data$gap) |>
    dplyr::mutate(pos = dplyr::row_number(),
                  pos = ifelse(.data$gap, NA, .data$pos)) |>
    dplyr::ungroup() |>
    dplyr::rename(!!paste0(col_name, '_pos') := 'pos') |>
    dplyr::select(-'gap')

}

