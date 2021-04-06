#' Multidimensional scaling plot of distances between methylation 
#' proportions (beta values)
#'
#' Same as \code{\link{plotMDS}}, except for an arc-sine transformation of the
#' methylation proportions.
#'
#' @param x \code{RangedSummarizedExperiment}, output from
#'   \code{calc_derivedasm} or \code{calc_asm}.
#' @param group Vector of group or any other labels, same length as number of
#'   samples.
#' @param top Number of top CpG sites used to calculate pairwise distances.
#' @param coverage Minimum number of reads covering a CpG site on each allele.
#'   Default = 5.
#' @param adj Text adjustment in y-axis. Default = 0.2.
#' @param pointSize Default = 4.
#'
#' @return Two-dimensional MDS plot.
#'
#' @importFrom SummarizedExperiment assays
#' @import ggplot2
#' 
#' @examples 
#' data(readtuples_output)
#' ASM <- calc_asm(readtuples_output)
#' grp <- factor(c(rep('CRC',3),rep('NORM',2)), levels = c('NORM', 'CRC'))
#' methyl_MDS_plot(ASM, grp)
#'
#' @export
methyl_MDS_plot <- function(x, group, top = 1000, coverage = 5, 
    adj = 0.02, pointSize = 4) {
    
    
    if (names(assays(x))[1] == "asm") {
        
        asm <- assays(x)[["asm"]]
        asm.red <- asm[rowSums(!is.na(assays(x)[["cov"]]) & 
            assays(x)[["cov"]] >= 
            coverage) == BiocGenerics::ncol(x), ]
        
        #mds_meth <- limma::plotMDS(asm.red, top = top, 
        #                            plot = FALSE)$cmdscale.out
        mds_meth <- limma::plotMDS(asm.red, top = top, 
                                    plot = FALSE)$eigen.vectors
        
    } else {
        
        full.cov <- assays(x)[["ref.cov"]] + assays(x)[["alt.cov"]]
        filt <- BiocGenerics::rowSums(full.cov >= coverage & 
            !is.na(full.cov)) == BiocGenerics::ncol(x)
        xfilt <- x[filt, ]
        
        # transform derived ASM
        methsTR <- as.matrix(assays(xfilt)[["der.ASM"]])
        
        methsTR <- asin(2 * methsTR - 1)
        bad <- BiocGenerics::rowSums(is.finite(methsTR)) < ncol(methsTR)
        if (any(bad)) 
            methsTR <- methsTR[!bad, , drop = FALSE]
        
        mds_meth <- limma::plotMDS(methsTR, top = top, 
                                    plot = FALSE)$eigen.vectors
    }
    
    df <- data.frame(dim1 = mds_meth[, 1], dim2 = mds_meth[, 
        2], names = colnames(x), treat = group)
    
    ggplot() + geom_point(data = df, mapping = aes_(x = ~dim1, 
        y = ~dim2, color = ~treat), size = pointSize) + geom_text(data = df, 
        mapping = aes_(x = ~dim1, y = ~dim2 - adj, label = ~names)) + 
        theme_bw()
    
}
