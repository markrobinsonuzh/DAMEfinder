#' Title
#'
#' @param sm_t
#' @param design
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_tstats <- function(sm_t, design, ...) {
  # moderated t-statistic
  coeff <- 2 # the column in the design matrix to consider
  fit <- limma::lmFit(sm_t, mod, method = "robust")
  fit2 <- limma::eBayes(fit, proportion = 0.01)
  fit2$t[, coeff]
}
