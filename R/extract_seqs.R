#' Extract founder (and reference) sequence from an alignment
#'
#' @param aln Alignment
#' @param founder_name Name of founder sequence in the alignment
#' @param ref_name Optional name of reference sequence in the alignment.
#' This can be used as input to the `ref_seq` argument in `run_wavess()`
#' (default: NULL, i.e. no reference sequence is returned)
#' @param start Start of founder sequence in alignment
#' (default: 1, i.e. beginning of sequence)
#' @param end End of founder sequence in alignment
#' (default: NULL, i.e. end of sequence)
#'
#' @return List of founder sequence and optional reference sequence
#' as character strings
#' @export
#'
#' @examples
#' hxb2_cons_founder <- ape::as.matrix.DNAbin(hxb2_cons_founder)
#' extract_seqs(hxb2_cons_founder, 'B.US.2011.DEMB11US006.KC473833', start = 6225, end = 7757)
#' extract_seqs(hxb2_cons_founder,
#' 'B.US.2011.DEMB11US006.KC473833', 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455',
#' start = 6225, end = 7757)
extract_seqs <- function(aln, founder_name, ref_name = NULL, start = 1, end = NULL){
  check_extract_seqs_inputs(aln, founder_name, ref_name, start, end)
  aln <- as.matrix(aln)
  if(is.null(end)){
    end <- ncol(aln)
  }
  ref <- NULL
  if(is.null(ref_name)){
    seqs <- founder_name
  }else{
    seqs <- c(founder_name, ref_name)
  }
  aln <- slice_aln(aln, start, end, seqs)
  aln <- aln[,as.character(aln[founder_name,,drop=TRUE]) != '-']
  founder <- toupper(paste0(as.character(aln[founder_name,]), collapse = ''))
  if(!is.null(ref_name)){
    ref <- toupper(paste0(as.character(aln[ref_name,]), collapse = ''))
  }
  return(list(founder=founder,ref=ref))
}

#' Slice alignment
#'
#' @param aln alignment to slice
#' @param start start position in alignment
#' @param end end position in alignment
#' @param seqs sequences to keep (default: labels(aln), i.e. all sequences)
#'
#' @return sliced alignment in `ape::DNAbin` format
#' @export
#'
#' @examples
#' slice_aln(hxb2_cons_founder, 1, 100)
slice_aln <- function(aln, start, end, seqs = labels(aln)){
  check_slice_aln_inputs(aln, start, end, seqs)
  return(as.matrix(aln)[seqs,start:end])
}

