#' Calculate region-level p-value
#'
#' This function uses the Simes method to calculate a regional-level p-value
#' based on the single-eBayes p-values. It goes into \code{\link{regionFinder}
#' to summarize each region detected.
#'
#' For positively correlated P-values, Simes method is even closer to the
#' nominal alpha level than the Bonferroni-Holm method.
#'
#' @param pval List of p-values obtained from \code{\link{get_tstats}}.
#' @return Vector of summarized pvals
#'
#' @export
simes_pval <- function(pval){
  #the simes thing for every cluster, pass on to regionFinder
  min((length(pval)*pval[order(pval)])/1:length(pval))
}

