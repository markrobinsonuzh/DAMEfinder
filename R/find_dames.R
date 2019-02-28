#' Find DAMEs
#'
#' This function finds Differential Allele-specific MEthylated regions (DAMEs).
#' It uses the \code{\link{regionFinder}} function from \code{bumphunter}.
#' 
#' A region is set as DAME if the smoothed t-Statistic (from \code{\link{eBayes}}) exceeds a cutoff K. 
#' K is set as the value of the 0.9 quantile of the absolute t-stats vector. When the smoothed t-stats 
#' in consecutive CpGs go above K or below -K, the region is set as a DAME.
#'
#' @param Q The percentile set to get a cutoff value K. K is the value on the Qth quantile
#' of the absolute values of the given smoothed t-Statistics vector. The default is set to 0.9.
#' @param maxGap Maximum gap between CpGs in a cluster (in bp). Default = 20.
#' @param verbose If the function should be verbose.
#' @param sa A SummarizedExperiment containing ASM values where each row and column correspond to a 
#' tuple/site and sample respectively.
#' @param control Index of columns corresponding to control samples.
#' @param treat Index of columns corresponding to treatment samples.
#' @param samp.perc At least this percentage of samples per group should have data 
#' (no NAs, coverage above threshold) per row. Default = 0.9.
#' @param coverage Minimum number of reads covering a CpG site (sum of both alleles) or tuple. 
#' Default = 5.
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
#' - pvalEmp: Empirical p-value obtained from permuting covariate of interest (group Vs treat).
#' @md
#' 
#'
#' @export
#'
#' @examples
find_dames <- function(sa, control, treat, samp.perc = 0.9, coverage = 5, Q=0.9, maxGap=20, verbose=TRUE, ...){
  
  pre_sa <- sa
  
  #get tstats with limma
  sat <- get_tstats(pre_sa, control, treat, 
                    samp.perc = samp.perc, 
                    coverage = coverage, 
                    maxGap = maxGap,
                    ...)

  # Detect real dames
  sm_tstat <- S4Vectors::mcols(sat)$smooth_tstat

  if(is.null(dim(SummarizedExperiment::assays(sat)[["asm"]]))){
    midpt <- BiocGenerics::start(sat)
  } else {
    midpt <- S4Vectors::mcols(sat)$midpt
  }

  K <- stats::quantile(abs(sm_tstat), Q, na.rm=TRUE)
  rf <- bumphunter::regionFinder(x = sm_tstat,
                                 chr = as.character(GenomeInfoDb::seqnames(sat)),
                                 pos = midpt,
                                 cluster = S4Vectors::mcols(sat)$cluster,
                                 cutoff = K,
                                 maxGap = maxGap,
                                 verbose = verbose)

  rf <- rf[,c("chr", "start", "end", "area")]
  rownames(rf) <- paste0("DAME", seq(from = 1, to = nrow(rf), by = 1))
  if(verbose) message(nrow(rf), " DAMEs found.")
  
  #### Permutation tests ####
  tot_length <- length(treat) + length(control)
  combs <- combn(tot_length, length(control))
  
  #Remove redundant perms
  if(length(control) == length(treat)){
  combs <- combs[,1:(ncol(combs)/2)]
  }
  
  #remove the true constrats
  rem <- apply(combs, 2, function(i){
    all.equal(i,control)
  })
  rem <- which(rem == T)
  if(!identical(rem, integer(0))) combs <- combs[,-rem]
  
  #Detect permuted dames
  message("Generating permutations")
  areas <- apply(combs, 2, function(i){
    treat_perm <- which(!(1:tot_length %in% i))
    sa_perm <- get_tstats(pre_sa, control = i, treat = treat_perm, 
                          samp.perc = samp.perc, 
                          coverage = coverage, 
                          maxGap = maxGap,
                          verbose = F,
                          ...)
    
    sm_tstat <- S4Vectors::mcols(sa_perm)$smooth_tstat
    
    if(is.null(dim(SummarizedExperiment::assays(sa_perm)[["asm"]]))){
      midpt <- BiocGenerics::start(sa_perm)
    } else {
      midpt <- S4Vectors::mcols(sa_perm)$midpt
    }
    
    K <- stats::quantile(abs(sm_tstat), Q, na.rm = TRUE)
    rf <- bumphunter::regionFinder(x = sm_tstat,
                                   chr = as.character(GenomeInfoDb::seqnames(sa_perm)),
                                   pos = midpt,
                                   cluster = S4Vectors::mcols(sa_perm)$cluster,
                                   cutoff = K,
                                   maxGap = maxGap,
                                   verbose = F)
    
    rf$area
  })
  
  all_areas <- sort(unlist(areas))
  total_areas <- length(all_areas)
  
  rf$pvalEmp <- sapply(rf$area, function(a){
    pperm <- (sum(all_areas > a) + 1) / (total_areas + 1)
  })
  
  return(rf)
}
