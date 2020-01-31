#' Find DAMEs
#'
#' This function finds Differential Allele-specific MEthylated regions (DAMEs).
#' It uses the \code{\link{regionFinder}} function from \code{bumphunter}, and
#' asigns p-values either empirically or using the Simes method.
#'
#' The simes method has higher power to detect DAMEs, but the consistency in
#' signal across a region is better controlled with the empirical method, since
#' it uses \code{regionFinder} and \code{getSegments} to find regions with
#' t-statistics above a cuttof (controled with parameter \code{Q}), whereas
#' with the 'simes' option, we initially detects clusters of CpG sites/tuples,
#' and then test if at least 1 differential site/tuple is present in the
#' cluster.
#'
#' We recommend trying out different \code{maxGap} and \code{Q} parameters,
#' since the size and the effect-size of obtained DAMEs change with these
#' parameters.
#'
#' @param sa A \code{SummarizedExperiment} containing ASM values where each row
#'   correspond to a tuple/site and a column to sample/replicate.
#' @param design A design matrix created with \code{\link{model.matrix}}.
#' @param coef Column in \code{design} specifying the parameter to estimate.
#'   Default = 2.
#' @param contrast a contrast matrix, generated with
#'   \code{\link{makeContrasts}}.
#' @param Q The percentile set to get a cutoff value K. K is the value on the
#'   Qth quantile of the absolute values of the given (smoothed) t-statistics.
#'   Only necessary if \code{pvalAssign} = 'empirical'. Default = 0.5.
#' @param pvalAssign Choose method to assign pvalues, either 'simes' (default)
#'   or 'empirical'. This second one performs \code{maxPerms} number of
#'   permutations to calculate null statistics, and runs \code{regionFinder}.
#' @param maxGap Maximum gap between CpGs in a cluster (in bp). NOTE: Regions
#' can be as small as 1 bp. Default = 20.
#' @param smooth Whether smoothing should be applied to the t-Statistics.
#'   Default = TRUE.
#' @param maxPerms Maximum possible permutations generated. Only necessary if
#' \code{pvalAssign} = 'empirical'. Default = 10.
#' @param method The method to be used in limma's \code{\link{lmFit}}. The
#'   default is set to 'ls' but can also be set to 'robust', which is
#'   recommended on a real data set.
#' @param trend Passed to \code{\link{eBayes}}. Should an intensity-trend be
#' allowed for the prior variance? Default is that the prior variance is
#' constant, e.g. FALSE.
#' @param verbose If the function should be verbose. Default = TRUE.
#' @param ... Arguments passed to \code{\link{get_tstats}}.
#'
#' @return A data frame of detected DAMEs ordered by the p-value. Each row
#'   is a DAME and the following information is provided in the columns
#'   (some column names change depending on the \code{pvalAssign} choice):
#'
#'   - chr: on which chromosome the DAME is found.
#'   - start: The start position of the DAME.
#'   - end: The end position of the DAME.
#'   - pvalSimes: p-value calculated with the Simes method.
#'   - pvalEmp: Empirical p-value obtained from permuting covariate of interest.
#'   - sumTstat: Sum of t-stats per segment/cluster.
#'   - meanTstat: Mean of t-stats per segment/cluster.
#'   - segmentL: Size of segmented cluster (from \code{\link{getSegments}}).
#'   - clusterL: Size of original cluster (from \code{\link{clusterMaker}}).
#'   - FDR: Adjusted p-value using the method of Benjamini, Hochberg. (from
#'   \code{\link{p.adjust}}).
#'   - numup: Number of sites with ASM increase in cluster (only for Simes).
#'   - numdown: Number of sites with ASM decrease in cluster (only for Simes).
#'
#' @md
#'
#' @examples
#' data(readtuples_output)
#' ASM <- calc_asm(readtuples_output)
#' grp <- factor(c(rep('CRC',3),rep('NORM',2)), levels = c('NORM', 'CRC'))
#' mod <- model.matrix(~grp)
#' dames <- find_dames(ASM, mod, verbose = FALSE)
#'
#' @export
#'
#'
find_dames <- function(sa, design, coef = 2, contrast = NULL, 
    smooth = TRUE, Q = 0.5, pvalAssign = "simes", maxGap = 20, 
    verbose = TRUE, maxPerms = 10, method = "ls", trend = FALSE, 
    ...) {
    
    
    if ("asm" %in% names(assays(sa))) 
        message("Using ASMtuple", appendLF = TRUE)
    if ("der.ASM" %in% names(assays(sa))) 
        message("Using ASMsnp", appendLF = TRUE)
    
    pre_sa <- sa
    
    # get tstats with limma
    sat <- get_tstats(pre_sa, design, maxGap = maxGap, coef = coef, 
        contrast = contrast, smooth = smooth, method = method, 
        trend = trend, ...)
    
    # choose smoothed if true
    if (smooth) {
        sm_tstat <- S4Vectors::mcols(sat)$smooth_tstat
    } else {
        sm_tstat <- S4Vectors::mcols(sat)$tstat
    }
    
    # choose position to find regions
    if (names(assays(sa))[1] == "asm") {
        midpt <- S4Vectors::mcols(sat)$midpt
    } else {
        midpt <- BiocGenerics::start(sat)
    }
    
    message("Detecting DAMEs", appendLF = TRUE)
    # detect dames
    K <- stats::quantile(abs(sm_tstat), Q, na.rm = TRUE)
    
    if (pvalAssign == "simes") {
        
        rf <- simes_pval(sat, sm_tstat, midpt)
        rf <- rf[order(rf$pvalSimes), ]
        rf$FDR <- stats::p.adjust(rf$pvalSimes, method = "BH")
        
    } else if (pvalAssign == "empirical") {
        
        rf <- bumphunter::regionFinder(x = sm_tstat, 
            chr = as.character(GenomeInfoDb::seqnames(sat)), 
            pos = midpt, cluster = S4Vectors::mcols(sat)$cluster, 
            cutoff = K, maxGap = maxGap, assumeSorted = TRUE, 
            order = FALSE, verbose = verbose)
        
        rf$pvalEmp <- empirical_pval(presa = pre_sa, design = design, 
            rforiginal = rf, coeff = coef, cont = contrast, smooth = smooth, 
            maxPerms = maxPerms, Q = Q, maxGap = maxGap, method = method, 
            ...)
        
        rf <- rf[order(rf$pvalEmp), ]
        
        rf$FDR <- stats::p.adjust(rf$pvalEmp, method = "BH")
        
        rf <- rf[, -c(6:8)]
        colnames(rf) <- c("chr", "start", "end", "meanTstat", 
            "sumTstat", "segmentL", "clusterL", "pvalEmp", "FDR")
        
    }
    if (verbose) 
        message(nrow(rf), " DAMEs found.", appendLF = TRUE)
    return(rf)
}
