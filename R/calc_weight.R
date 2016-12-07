# calculate the weight per site given beta and a
#' Title
#'
#' @param beta
#' @param MM
#' @param UU
#' @param a
#'
#' @return
#' @export
#'
#' @examples
#'
calc_weight <- function(beta, MM, UU, a) {

  # set quantile range
  theta = seq(0.005, 0.995, length = 1000)

  # calculate pbeta (note that this is the cumulative prob distr). To get the densities, use dbeta.
  posterior <- pbeta(theta, beta+MM, beta+UU, lower.tail = T)

  # get prob(0.3<b<0.7, ie a=0.2): our new weight
  w_1 <- which(theta<(0.5-a))
  index_1 <- w_1[length(w_1)] + 1

  w_2 <- which(theta<(0.5+a))
  index_2 <- w_2[length(w_2)]

  posterior[index_2] - posterior[index_1]

}


