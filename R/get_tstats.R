#' Get t-Statistics
#'
#' This function calculates a moderated t-Statistic per site or tuple using limma's lmFit and eBayes
#' functions. It then smoothes the obtained t-Statistics using bumphunter's smoother function.
#' The smoothing is done on genomic clusters consisting of CpGs that are close to each other.
#' In the case of tuples, the midpoint of the two genomic positions in each tuple is used as the
#' genomic position of that tuple, to apply smoothing on consecutive positions (or tuples).
#' The function takes a matrix containing ASM across samples and a design matrix in which the
#' sample conditions are specified. The t-Statistic reflects a measure of difference in
#' ASM scores between the different sample conditions.
#'
#' @param sa A SummarizedExperiment containing ASM values where each row and column correspond to a tuple/site and
#' sample respectively.
#' @param control Index of columns corresponding to control samples.
#' @param treat Index of columns corresponding to treatment samples.
#' @param samp.perc At least this percetage of samples per group should have data (no NAs, coverage above threshold)
#' per row. Default = 0.9.
#' @param coverage Minimum number of reads covering a CpG site/tuple. Default = 5.
#' @param method The method to be used in limma's lmFit. The default is set to "ls" but one can
#' also set it to "robust", which is recommended on a real data set.
#' @param maxGap The maximum allowed gap between genomic positions for clustering of genomic
#' regions to be used in smoothing. Default = 20.
#' @param coef Columm taken from limma "fit".
#' @param verbose Set verbose. Default = TRUE.
#'
#' @return A vector of smoothed t-Statistics within the SummarizedExperiment
#' @importFrom BiocGenerics start
#' @importFrom SummarizedExperiment assays
#'
#' @export
#' @examples

get_tstats <- function(sa, control, treat, samp.perc = 0.9, coverage = 5, method="ls", maxGap=20, coef=2, verbose=TRUE) {

  # choose SumExp type and filter
  if(dim(SummarizedExperiment::rowData(sa))[2] > 0){

    control.covfilt <- BiocGenerics::rowSums(assays(sa)[["cov"]][,control] >= coverage &
                                 !is.na(assays(sa)[["cov"]][,control])) >= floor(length(control)*samp.perc)

    treat.covfilt <- BiocGenerics::rowSums(assays(sa)[["cov"]][,treat] >= coverage &
                               !is.na(assays(sa)[["cov"]][,treat])) >= floor(length(treat)*samp.perc)

    sa <- sa[control.covfilt&treat.covfilt,]
    asm <- assays(sa)[["asm"]]

  } else {

    full.cov <- assays(sa)[["ref.cov"]] + assays(sa)[["alt.cov"]]

    control.covfilt <- BiocGenerics::rowSums(full.cov[,control] >= coverage &
                                               !is.na(full.cov[,control])) >= floor(length(control)*samp.perc)

    treat.covfilt <- BiocGenerics::rowSums(full.cov[,treat] >= coverage &
                                             !is.na(full.cov[,treat])) >= floor(length(treat)*samp.perc)

    sa <- sa[control.covfilt&treat.covfilt,]
    asm <- assays(sa)[["der.ASM"]]
  }

  # set design matrix
  cols <- c(rep(1,length(treat)),rep(0,length(control)))
  mod <- matrix(data = c(rep(1,BiocGenerics::ncol(sa)), cols), ncol = 2)

  # moderated t-statistic using specified column in the design matrix
  if(verbose) message("Calculating moderated t-statistics")
  fit <- limma::lmFit(asm, mod, method = method)
  fit2 <- limma::eBayes(fit)
  S4Vectors::mcols(sa)$tstat <- fit2$t[, coef]

  # smooth moderated t-stats
  if(verbose) message("Smoothing moderated t-statistics")

  if(dim(SummarizedExperiment::rowData(sa))[2] > 1){
    midpt <- S4Vectors::mcols(sa)$midpt
  } else {
    midpt <- start(sa)
  }

  S4Vectors::mcols(sa)$cluster <- pns <- bumphunter::clusterMaker(as.character(GenomeInfoDb::seqnames(sa)),
                                                                  midpt, maxGap=maxGap)

  smooth <- bumphunter::smoother(y = S4Vectors::mcols(sa)$tstat, x = midpt,
                                cluster = pns,
                                smoothFunction = bumphunter::loessByCluster,
                                verbose = verbose)

  S4Vectors::mcols(sa)$smooth_tstat <- smooth$fitted[,1]
  sa
}
