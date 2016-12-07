#' Title
#'
#' @param nm
#' @param smoothed_betas
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
find_dames <- function(nm, smoothed_betas, verbose=TRUE, ...) {

  chr <- limma::strsplit2(nm, ".", fixed=T)[,1]
  pos1 <- as.numeric(limma::strsplit2(nm, ".", fixed=T)[,2])
  pos2 <- as.numeric(limma::strsplit2(nm, ".", fixed=T)[,3])
  midpt <- floor((pos2 - pos1)/2)
  pos <- pos1 + midpt

  pns <- bumphunter::clusterMaker(chr, pos, maxGap = 300)

  Q <- 0.9
  # Detect DAMEs
  K <- quantile(abs(smoothed_betas),Q, na.rm=T)
  dames <- bumphunter::regionFinder(x = smoothed_betas, chr = chr, pos = pos, cluster = pns,
                        cutoff = K, verbose = verbose)


}







