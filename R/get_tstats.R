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
#' @param betas
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_smoothed_tstats <- function(betas, Q=0.9, verbose=TRUE, maxGap=300) {
  
  rows <- names(betas)
  chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
  pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
  pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])
  midpt <- floor((pos2 - pos1)/2)
  pos <- pos1 + midpt
  
  pns <- bumphunter::clusterMaker(chr, pos, )
  
  smooth <- bumphunter::smoother(y = betas, x = pos, cluster = pns, 
                                 smoothFunction = loessByCluster,
                                 verbose = verbose)
  smooth$fitted
}
