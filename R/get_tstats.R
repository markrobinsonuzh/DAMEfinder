#' Title
#'
#' @param method 
#' @param coef 
#' @param design
#' @param sa 
#' @param maxGap 
#'
#' @return
#' @export
#'
#' @examples
get_tstats <- function(sa, design, method="robust", maxGap=300, coef=2, verbose=TRUE) {
  
  asm <- SummarizedExperiment::assays(sa)[["asm"]]
  
  # moderated t-statistic using specified column in the design matrix
  if(verbose) message("Calculating moderated t-statistics.")
  fit <- limma::lmFit(asm, design, method = method)
  fit2 <- limma::eBayes(fit)
  S4Vectors::mcols(sa)$tstat <- fit2$t[, coef]

  # smooth moderated t-stats
  if(verbose) message("Smoothing moderated t-statistics.")
  midpt <- S4Vectors::mcols(sa)$midpt  
  S4Vectors::mcols(sa)$cluster <- pns <- bumphunter::clusterMaker(GenomeInfoDb::seqnames(sa), 
                                                                  midpt, maxGap=maxGap)
  
  smooth <- bumphunter::smoother(y = S4Vectors::mcols(sa)$tstat, x = midpt, 
                                 cluster = pns, 
                                 smoothFunction = bumphunter::loessByCluster,
                                 verbose = verbose)
  S4Vectors::mcols(sa)$smooth_tstat <- smooth$fitted[,1]
  sa
}
