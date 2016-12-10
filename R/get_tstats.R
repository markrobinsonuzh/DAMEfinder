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
  
  asm <- assays(sa)[["asm"]]
  
  # moderated t-statistic using specified column in the design matrix
  if(verbose) message("Calculating moderated t-statistics.")
  fit <- limma::lmFit(asm, design, method = method)
  fit2 <- limma::eBayes(fit)
  mcols(sa)$tstat <- fit2$t[, coef]

  # smooth moderated t-stats
  if(verbose) message("Smoothing moderated t-statistics.")
  midpt <- mcols(sa)$midpt  
  mcols(sa)$cluster <- pns <- bumphunter::clusterMaker(seqnames(sa), 
                                                       midpt, maxGap=maxGap)
  
  smooth <- bumphunter::smoother(y = mcols(sa)$tstat, x = midpt, 
                                 cluster = pns, 
                                 smoothFunction = loessByCluster,
                                 verbose = verbose)
  mcols(sa)$smooth_tstat <- smooth$fitted[,1]
  sa
}
