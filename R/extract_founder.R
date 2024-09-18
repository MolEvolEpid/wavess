#' Extract founder sequence from an alignment
#'
#' @param aln Alignment
#' @param founder_name Name of founder sequence in the alignment
#' @param start Start of founder sequence (default: 1, i.e. beginning of sequence)
#' @param end End of founder sequence (default: NULL, i.e. end of sequence)
#'
#' @return Founder sequence as a character string
#' @export
#'
#' @examples
#' hxb2_founder <- ape::as.matrix.DNAbin(hxb2_founder)
#' extract_founder(hxb2_founder, 'B.US.2011.DEMB11US006.KC473833')
extract_founder <- function(aln, founder_name, start = 1, end = NULL){
  check_extract_founder_inputs(aln, founder_name, start, end)
  aln <- as.matrix(aln)
  if(is.null(end)){
    end <- ncol(aln)
  }
  return(toupper(gsub('-', '', paste0(as.character(aln[founder_name,start:end,drop=TRUE]), collapse = ''))))
}

