#' Find DAMEs
#'
#' This function finds Differentially Allele-specifically MEthylated regions (DAMEs).
#' It uses the 'regionFinder' function from the bumphunter package. A region is set as
#' a DAME if the smoothed t-Statistic, which is set as input to the function, exceeds
#' a cutoff K. K is set as the value of the 0.9 quantile of the absolute t-Statistics vector.
#' When the smoothed t-Stats in consecutive CpGs go above K or below -K, the region
#' is set as a DAME.
#'
#' @param Q The percentile set to get a cutoff value K. K is the value on the Qth quantile
#' of the absolute values of the given smoothed t-Statistics vector. The default is set to 0.9.
#' @param maxGap Maximum gap between CpGs in a cluster.
#' @param verbose If the function should be verbose.
#' @param sa Vector of smoothed t-Statistics obtained from the get_stats function.
#'
#' @return A data frame of detected DAMEs ordered by the area each DAME has above the cutoff K.
#' The larger the area value, the more important the DAME. Each row refers to a DAME and the
#' following information is provided in the columns:
#'
#' - chr: on which chromosome the DAME is found
#' - start: The start position of the DAME
#' - end: The end position of the DAME
#' - value:
#' - area: The area of the DAME beyond the exceeded cutoff K.
#' - cluster: Then genomic cluster the DAME belongs to. Smoothing of the t-Statistics was
#' done per cluster.
#' - indexStart:
#' - indexEnd:
#' - L:
#' - clusterL:
#' @md
#'
#' @export
#'
#' @examples
find_dames <- function(sa, Q=0.9, maxGap=20, verbose=TRUE) {

  # Detect DAMEs
  sm_tstat <- S4Vectors::mcols(sa)$smooth_tstat
  #sm_tstat <- S4Vectors::mcols(sa)$tstat

  if(is.null(dim(SummarizedExperiment::assays(sa)[["asm"]]))){
    midpt <- BiocGenerics::start(sa)
  } else {
    midpt <- S4Vectors::mcols(sa)$midpt
  }

  K <- stats::quantile(abs(sm_tstat), Q, na.rm=TRUE)
  rf <- bumphunter::regionFinder(x = sm_tstat,
                                 chr = as.character(GenomeInfoDb::seqnames(sa)),
                                 pos = midpt,
                                 cluster = S4Vectors::mcols(sa)$cluster,
                                 cutoff = K,
                                 verbose = verbose)

  rf <- rf[,c("chr", "start", "end", "value", "area")]
  rownames(rf) <- paste0("DAME", seq(from = 1, to = nrow(rf), by = 1))
  if(verbose) message(nrow(rf), " DAMEs found.")
  rf
}







