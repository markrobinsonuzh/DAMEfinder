#' Find DAMEs
#'
#' This function finds Differential Allele-specific MEthylated regions (DAMEs).
#' It uses the \code{\link{regionFinder}} function from \code{bumphunter}.
#' 
#' A region is set as DAME if the smoothed t-Statistic (from \code{\link{eBayes}}) exceeds a cutoff K. 
#' K is set as the value of the 0.9 quantile of the absolute t-stats vector. When the smoothed t-stats 
#' in consecutive CpGs go above K or below -K, the region is set as a DAME.
#'
#'
#' @param sa A SummarizedExperiment containing ASM values where each row and column correspond to a 
#' tuple/site and sample respectively.
#' @param design A design matrix created with \code{\link{model.matrix}}.
#' @param coef Column in model.matrix specifying the parameter to estimate. Default = 2.
#' @param Q The percentile set to get a cutoff value K. K is the value on the Qth quantile
#' of the absolute values of the given smoothed t-Statistics vector. The default is set to 0.9.
#' @param pvalAssign Choose method to assign pvalues, either "simes" (default) or "empirical",
#' which performs a number of permutations to calculate the null.
#' @param maxGap Maximum gap between CpGs in a cluster (in bp). Default = 20.
#' @param smooth Whether smoothing should be applied to the t-Statistics. Default = TRUE
#' @param verbose If the function should be verbose.
#' @param maxPerms = 10
#' @param ... arguments passed to \code{\link{get_tstats}}.
#'
#' @return A data frame of detected DAMEs ordered by the area each DAME has above the cutoff K.
#' The larger the area value, the more important the DAME. Each row refers to a DAME and the
#' following information is provided in the columns:
#'
#' - chr: on which chromosome the DAME is found.
#' - start: The start position of the DAME.
#' - end: The end position of the DAME.
#' - area: Sum (abs) of all t-stat values in a region.
#' - pvalEmp: Empirical p-value obtained from permuting covariate of interest, 
#' or through Simes correction.
#'
#' @md
#' 
#' @export
#'
#' @examples
find_dames <- function(sa, design, coef = 2, smooth = TRUE, Q = 0.9, pvalAssign = "simes", 
                       maxGap = 300, verbose = TRUE, maxPerms = 10, ...){
  
  pre_sa <- sa
  
  #get tstats with limma
  sat <- get_tstats(pre_sa, design,
                    maxGap = maxGap,
                    coef = coef,
                    smooth = smooth,
                    ...)

  # choose smoothed if true
  if(smooth){
  sm_tstat <- S4Vectors::mcols(sat)$smooth_tstat
  } else {
    sm_tstat <- S4Vectors::mcols(sat)$tstat
  }

  #choose position to find regions
  if(names(assays(sa))[1] == "asm"){
    midpt <- S4Vectors::mcols(sat)$midpt
  } else {
    midpt <- BiocGenerics::start(sat)
  }

  message("Detecting DAMEs", appendLF = TRUE)
  #detect dames
  K <- stats::quantile(abs(sm_tstat), Q, na.rm=TRUE)
  
  if(pvalAssign == "simes"){
    
    rf <- bumphunter::regionFinder(x = sm_tstat,
                                   chr = as.character(GenomeInfoDb::seqnames(sat)),
                                   pos = midpt,
                                   cluster = S4Vectors::mcols(sat)$cluster, #This might make no difference?
                                   y = S4Vectors::mcols(sat)$p.value,
                                   summary = simes_pval,
                                   cutoff = K,
                                   maxGap = maxGap,
                                   assumeSorted = TRUE,
                                   order = FALSE,
                                   verbose = verbose)
    
    rf <- rf[,-c(6:8)]
    colnames(rf) <- c("chr","start","end","pvalSimes","sumPvalues","segmentL","clusterL")
    rf <- rf[order(rf$pvalSimes),]
    rf$FDR <- p.adjust(rf$pvalSimes, method="BH")
    
  } else if(pvalAssign == "empirical"){
    
    rf <- bumphunter::regionFinder(x = sm_tstat,
                                   chr = as.character(GenomeInfoDb::seqnames(sat)),
                                   pos = midpt,
                                   cluster = S4Vectors::mcols(sat)$cluster, #This might make no difference?
                                   cutoff = K,
                                   maxGap = maxGap,
                                   assumeSorted = TRUE,
                                   verbose = verbose)
    
    rf <- empirical_pval(pre_sa = pre_sa, 
                         design = design, 
                         rforiginal = rf, 
                         coeff = coef, 
                         smooth = smooth, 
                         maxPerms = maxPerms, 
                         K = K, 
                         maxGap = maxGap, 
                         ...)
    
  }
  if(verbose) message(nrow(rf), " DAMEs found.")
  return(rf)
}
