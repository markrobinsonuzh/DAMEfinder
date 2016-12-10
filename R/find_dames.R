#' Title
#'
#' @param Q 
#' @param maxGap 
#' @param verbose
#' @param sa 
#'
#' @return
#' @export
#'
#' @examples
find_dames <- function(sa, Q=0.9, maxGap=300, verbose=TRUE) {

  # Detect DAMEs
  sm_tstat <- mcols(sa)$smooth_tstat
  K <- quantile(sm_tstat, Q, na.rm=TRUE)
  rf <- bumphunter::regionFinder(x = sm_tstat, chr = as.character(GenomeInfoDb::seqnames(sa)), 
                           pos = mcols(sa)$midpt, cluster = mcols(sa)$cluster,
                           cutoff = K, verbose = verbose)
  if(verbose) message(nrow(rf), " DAMEs found.")
  rf
}







