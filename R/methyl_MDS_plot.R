#' Multidimensional scaling plot of distances between methylation proportions
#'
#' Same as limma::plotMDS, expect for an arc-sine transformation of the methylation proportions
#' or scores from 0 to 1.
#' TODO: check what transformation for ASM-score
#'
#' @param x SummarizedExperiment, output from calc_derivedasm(). TODO: also for calc_asm()?
#' @param top Number of top CpG sites used to calculate pairwise distances.
#' @param coverage Minimum number of reads covering a CpG site on each allele. Default = 4.
#' @param color Vector of group labels same length as number of samples
#' #TRUE or FALSE. the object should have group informationin colData
#'
#' @return Plot
#' @examples
#'
#'
#' @importFrom SummarizedExperiment assays
#' @import ggplot2
#'
#' @export
#'
#'
#'
methyl_MDS_plot <- function(x, top = 1000, coverage = 4, color){

  if(dim(SummarizedExperiment::rowData(x))[2] > 0){

    asm <- assays(x)[["asm"]]
    asm.red <- asm[rowSums(!is.na(assays(x)[["cov"]]) & assays(x)[["cov"]] >= coverage) == BiocGenerics::ncol(x),]
    mds_meth <- limma::plotMDS(asm.red)$cmdscale.out

  } else {
    #filter ref and alt coverage
    # ref.cov <- BiocGenerics::rowSums(assays(x)[["ref.cov"]] >= coverage &
    #                              !is.na(assays(x)[["ref.cov"]])) >= BiocGenerics::ncol(x)
    #
    # alt.cov <- BiocGenerics::rowSums(assays(x)[["alt.cov"]] >= coverage &
    #                            !is.na(assays(x)[["alt.cov"]])) >= BiocGenerics::ncol(x)
    #
    # xfilt <- x[ref.cov&alt.cov,]

    full.cov <- assays(x)[["ref.cov"]] + assays(x)[["alt.cov"]]
    filt <- BiocGenerics::rowSums(full.cov >= coverage & !is.na(full.cov)) == BiocGenerics::ncol(x)
    xfilt <- x[filt,]

    #transform derived ASM
    methsTR <- as.matrix(assays(xfilt)[["der.ASM"]])

    methsTR <- asin(2*methsTR-1)
    bad <- BiocGenerics::rowSums(is.finite(methsTR)) < ncol(methsTR)
    if(any(bad)) methsTR <- methsTR[!bad,,drop=FALSE]

    #Grab top sites
    nprobes <- BiocGenerics::nrow(methsTR)
    nsamples <- BiocGenerics::ncol(methsTR)
    top <- min(top,nprobes)

    #Distance measure is mean of top squared deviations for each pair of arrays
    topindex <- nprobes-top+1L
    dd <- matrix(0,nrow=nsamples,ncol=nsamples)
    for (i in 2L:(nsamples))
      for (j in 1L:(i-1L))
        dd[i,j] <- sqrt(mean(sort.int((methsTR[,i]-methsTR[,j])^2,partial=topindex)[topindex:nprobes]))

    #multidimensional scaling
    mds_meth <- stats::cmdscale(stats::as.dist(dd),2)
  }

  df <- data.frame(dim1 = mds_meth[,1], dim2 = mds_meth[,2], names = colnames(x), treat = color)

  ggplot(df)+
    geom_point(aes(x=dim1, y=dim2, color = treat), size=5) +
    geom_text(aes(x=dim1, y=dim2-.02, label = names), size=4)+
    theme_bw()
}
