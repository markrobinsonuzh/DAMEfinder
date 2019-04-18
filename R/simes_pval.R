#' Calculate region-level p-value
#'
#' This function uses the Simes method to calculate a regional-level p-value
#' based on the single-eBayes p-values. It highly depends on the choice of
#' \code{maxGap} in \code{find_dames}.
#'
#' When used as a FDR-control method, for positively correlated P-values, Simes
#' method is even closer to the nominal alpha level than the Bonferroni-Holm
#' method.
#'
#' @param sat Output from \code{\link{get_tstats}}.
#' @param smtstat (Smoothed) tstat vector from \code{\link{get_tstats}}.
#' @param midpt Coordinate vector for each CpG site/tuple.
#' @return Vector of summarized pvals
#'
#' @export
simes_pval <- function(sat, smtstat, midpt){

  #try to use this strategy using full clusters, not segments which means
  #choosing the maxGap for cluster in get_tstat is the most important parameter
  #for pval method
  
  pval <- S4Vectors::mcols(sat)$p.value
  cluster <- mcols(sat)$cluster
  cluster.ids <- unique(cluster)
  
  chr <- as.character(seqnames(sat))
  starts <- midpt
  ends <- midpt
  
  simes <- function(pval) min((length(pval)*pval[order(pval)])/1:length(pval))

  realregs <- data.frame(
    chr=sapply(cluster.ids,
               function(Index) chr[cluster == Index][1]),
    start=sapply(cluster.ids,
                 function(Index) min(starts[cluster == Index])),
    end=sapply(cluster.ids, 
               function(Index) max(ends[cluster == Index])),
    meanTstat=sapply(cluster.ids, 
                     function(Index) mean(smtstat[cluster == Index])),
    sumTstat=sapply(cluster.ids, 
                    function(Index) sum(smtstat[cluster == Index])),
    pvalSimes=sapply(cluster.ids, 
                     function(Index) simes(pval[cluster == Index])),
    clusterL=sapply(cluster.ids, 
                 function(Index) length(cluster[cluster == Index])),
    numup=sapply(cluster.ids, 
                 function(Index) sum(smtstat[cluster == Index] > 0)),
    numdown=sapply(cluster.ids, 
                   function(Index) sum(smtstat[cluster == Index] < 0)))
  
  return(realregs)
}