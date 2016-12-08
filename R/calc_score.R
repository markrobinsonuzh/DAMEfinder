
#' Get a vector of the posterior scores for all positions given beta and a
#'
#' @param beta
#' @param a
#' @param df
#'
#' @return
#' @export
#'
#' @examples
calc_score <- function(beta, a, df) {

  weights <- parallel::mcmapply(calcWeight, beta, df$MM, df$UU, a, mc.cores = min(7, ceiling(parallel::detectCores()/6)))
  df$log_odds_ratio*weights
}

#' Calculate the odds ratio
#'
#' @param s
#'
#' @return
#' @export
#'
#' @examples
calc_logodds <- function(s, eps=1) {
  ratio <- with(s, ( (MM+eps)*(UU+eps) ) / ( (MU+eps)*(UM+eps) ) )
  s$logodds <- log10(ratio)
  s
}



#' Transform the ASM scores
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
transform_scores <- function(x) {
  modulus_sqrt(x)
}



#' Function to calculate signed square root (aka modulus square root)
#'
#' @param values
#'
#' @return
#' @export
#'
#' @examples
modulus_sqrt <- function(values) {
  sign(values)*sqrt(abs(values))
}



