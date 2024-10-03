#' Calculate diversity and divergence
#'
#' See `[vignette('analyze_output')]` for more details.
#'
#' @param aln DNA alignment in `[ape::DNAbin]` format
#' @param founder Name of the founder sequence in the alignment
#' @param gen Vector that indicates the generation of each sequence in the alignment,
#' assumed to be in the same order as the alignment
#'
#' @return tibble including mean divergence for each generation
#' @export
#'
#' @examples
#' # This example is somewhat contrived, but it shows how it works.
#' calc_div_metrics(hxb2_cons_founder, 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', c(1,2,2))
calc_div_metrics <- function(aln, founder, gen){
  check_name_in_alignment(aln, founder, 'aln', 'founder')
  if(length(labels(aln)) != length(gen)){
    stop('The length of gen must be the same as the number of sequences in the alignment')
  }
  aln <- as.matrix(aln)
  distmat <- ape::dist.dna(aln, model = 'raw', as.matrix = TRUE)
  lapply(unique(gen), function(x){
    distmat_sub <- distmat[rownames(distmat) == founder | gen == x,
                           colnames(distmat) == founder | gen == x, drop = FALSE]
    divergence <- sum(distmat_sub[founder,])/sum(gen == x)
    distmat_sub <- distmat_sub[rownames(distmat_sub) != founder,
                           colnames(distmat_sub) != founder, drop = FALSE]
    diversity <- sum(distmat_sub[lower.tri(distmat_sub)])/choose(nrow(distmat_sub), 2)
    return(tibble::tibble(gen = x,
                          diversity = diversity,
                          divergence = divergence))
  }) |> dplyr::bind_rows()
}
