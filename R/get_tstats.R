#' Get t-Statistics on Tuples
#'
#' This function calculates a moderated t-Statistic per tuple using limma's lmFit and eBayes 
#' functions. It then smoothes the obtained t-Statistics using bumphunter's smoother function. 
#' The smoothing is done on genomic clusters consisting of CpGs that are close to each other.
#' The midpoint of the two genomic positions in each tuple is used as the genomic position of
#' that tuple, to apply smoothing on consecutive positions (or tuples). The function takes a
#' matrix containing ASM scores per tuple across samples and a design matrix in which the 
#' sample conditions are specified. The t-Statistic reflects a measure of difference in 
#' ASM scores between the different sample conditions for every tuple.
#'
#' @param method The method to be used in limma's lmFit. The default is set to "ls" but one can 
#' also set it to "robust", which is recommended on a real data set. 
#' @param coef 
#' @param design
#' @param sa The matrix containing ASM scores where each row and column correspond to a tuple and
#' sample respectively.
#' @param maxGap The maximum allowed gap between genomic positions for clustering of genomic
#' regions to be used in smoothing. The default is set to 300.
#'
#' @return A vector of smoothed t-Statistics.
#' @export
#'
#' @examples
#' tStatistics <- get_tstats(ASM_score_matrix, dsgn)
get_tstats <- function(sa, design, method="ls", maxGap=300, coef=2, verbose=TRUE) {
  
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
