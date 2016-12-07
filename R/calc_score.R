
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
process <- function(s) {
  s$strand <- "*"
  s$MM <- s$MM + 1
  s$MU <- s$MU + 1
  s$UM <- s$UM + 1
  s$UU <- s$UU + 1
  s$cov <- s$cov + 4
  # calc log odds ratio
  ratio <- (s$MM*s$UU)/(s$MU*s$UM)
  s$log_odds_ratio <-  log(ratio, base=10)
  return(s)
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



