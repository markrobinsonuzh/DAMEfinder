#' Title
#'
#' @param betas
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_smoothed_t_stats <- function(betas, ...) {

  rows <- names(betas)
  chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
  pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
  pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])
  midpt <- floor((pos2 - pos1)/2)
  pos <- pos1 + midpt

  pns <- bumphunter::clusterMaker(chr, pos, maxGap = 300)

  verbose <- TRUE
  Q <- 0.9

  smooth <- bumphunter::smoother(y = betas, x = pos, cluster = pns, smoothFunction = loessByCluster,
                     verbose = verbose)
  smooth$fitted
}
