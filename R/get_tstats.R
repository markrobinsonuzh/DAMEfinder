#' Title
#'
#' @param sm_t
#' @param method 
#' @param coef 
#' @param design
#'
#' @return
#' @export
#'
#' @examples
get_tstats <- function(sm_t, design, method="robust", coef=2) {
  # moderated t-statistic
   # the column in the design matrix to consider
  fit <- limma::lmFit(sm_t, design, method = method)
  fit2 <- limma::eBayes(fit)
  fit2$t[, coef]
}


#' Title
#'
#' @param Q 
#' @param verbose 
#' @param maxGap 
#' @param betas
#'
#' @return
#' @export
#'
#' @examples
get_smoothed_tstats <- function(betas, Q=0.9, verbose=TRUE, maxGap=300) {
  
  rows <- names(betas)
  ss <- limma::strsplit2(rows, ".", fixed=TRUE)
  chr <- ss[,1]
  midpt <- (as.numeric(ss[,2])+as.numeric(ss[,3]))/2

  pns <- bumphunter::clusterMaker(chr, midpt, maxGap=maxGap)
  
  smooth <- bumphunter::smoother(y = betas, x = midpt, cluster = pns, 
                                 smoothFunction = loessByCluster,
                                 verbose = verbose)
  smooth$fitted
}
