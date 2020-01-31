#' Plot score tracks
#'
#'
#' @param dame GRanges object containing a region of interest, or detected with
#'   find_dames
#' @param window Number of CpG sites outside (up or down-stream) of the DAME
#'   should be plotted. Default = 0.
#' @param positions Number of bp sites outside (up or down-stream) of the DAME
#'   should be plotted. Default = 0.
#' @param derASM SummarizedExperiment object obtained from calc_derivedasm
#'   (Filtering should be done by the user)
#' @param ASM SummarizedExperiment object obtained from calc_asm (Filtering
#'   should be done by the user)
#' @param colvec Vector of colors (mainly useful for the SNP plot, because I add
#'   it with cowplot, so I don't export a ggplot, optional)
#' @param plotSNP whether to add the SNP track, only if derASM is specified.
#'   Default = FALSE
#'
#' @return Plot
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics start<-
#' @importFrom BiocGenerics end<-
#' @import ggplot2
#' 
#' @examples 
#' library(GenomicRanges)
#' DAME <- GRanges(19, IRanges(306443,310272))
#' data('readtuples_output')
#' ASM <- calc_asm(readtuples_output)
#' SummarizedExperiment::colData(ASM)$group <- c(rep('CRC',3),rep('NORM',2))
#' SummarizedExperiment::colData(ASM)$samples <- colnames(ASM)
#' dame_track(dame = DAME, 
#'            ASM = ASM) 
#'
#' @export

dame_track <- function(dame, window = 0, positions = 0, derASM = NULL, 
    ASM = NULL, colvec = NULL, plotSNP = FALSE) {
    
    res_dame <- dame
    start(res_dame) <- start(dame) - positions
    end(res_dame) <- end(dame) + positions
    
    if (!is.null(derASM)) {
        
        if (!all.equal(colData(derASM)$samples, colnames(derASM))) {
            stop("Sample names in colData() and colnames are different")
        }
        
        ASMsnp <- assay(derASM, "der.ASM")
        SNP <- assay(derASM, "snp.table")
        ref <- assay(derASM, "ref.meth")/assay(derASM, "ref.cov")
        alt <- assay(derASM, "alt.meth")/assay(derASM, "alt.cov")
        
        snpgr <- SummarizedExperiment::rowRanges(derASM)
        
        over <- GenomicRanges::findOverlaps(snpgr, res_dame)
        if (window != 0) {
            win <- c(seq(from = (queryHits(over)[1] - window), 
                to = (queryHits(over)[1] - 1), by = 1), queryHits(over), 
                seq(from = (utils::tail(queryHits(over), n = 1) + 
                    1), to = (utils::tail(queryHits(over), n = 1) + 
                    window), by = 1))
        } else {
            win <- queryHits(over)
        }
        
        # ASMsnp
        subASMsnp <- as.data.frame(ASMsnp[win, ])
        subref <- as.data.frame(ref[win, ])
        subalt <- as.data.frame(alt[win, ])
        
        subASMsnp$pos <- subref$pos <- subalt$pos <- start(snpgr)[win]
        subASMsnp_long <- reshape2::melt(subASMsnp, id.vars = "pos", 
            measure.vars = colnames(ASMsnp))
        subASMsnp_long$score <- "ASMsnp"
        
        # marg.meth per allele
        subref_long <- reshape2::melt(subref, id.vars = "pos", 
            measure.vars = colnames(ASMsnp))
        subref_long$score <- "REF:meth"
        
        subalt_long <- reshape2::melt(subalt, id.vars = "pos", 
            measure.vars = colnames(ASMsnp))
        subalt_long$score <- "ALT:meth"
        
        # SNP
        subSNP <- SNP[win, ]
        subSNP[is.na(subSNP)] <- "none"
        subSNP <- data.frame(subSNP, stringsAsFactors = FALSE)
        subSNP$pos <- start(snpgr)[win]
        subSNP_long <- reshape2::melt(subSNP, id.vars = "pos", 
            measure.vars = colnames(ASMsnp), value.name = "snp.pos")
        
        # Make SNP name nicer
        loc <- as.integer(stringr::str_extract(subSNP_long$snp.pos, 
            "[0-9]+$"))
        chrom <- as.integer(stringr::str_extract(subSNP_long$snp.pos, 
            "^[0-9]+"))
        subSNP_long$snp.pos <- ifelse(subSNP_long$snp.pos == 
            "none", "none", sprintf("chr%s:%s", chrom, loc))
    }
    
    if (!is.null(ASM)) {
        
        if (!all.equal(colData(ASM)$samples, colnames(ASM))) {
            stop("Sample names in colData() and colnames are different")
        }
        
        meth <- (assay(ASM, "MM") + assay(ASM, "MU") + assay(ASM, 
            "UM"))/assay(ASM, "cov")
        
        ASMtuple <- assay(ASM, "asm")
        tupgr <- SummarizedExperiment::rowRanges(ASM)
        
        
        # ASMtuple
        over <- GenomicRanges::findOverlaps(tupgr, res_dame)
        
        if (window != 0) {
            win <- c(seq(from = (queryHits(over)[1] - window), 
                to = (queryHits(over)[1] - 1), by = 1), queryHits(over), 
                seq(from = (utils::tail(queryHits(over), n = 1) + 
                1), to = (utils::tail(queryHits(over), n = 1) + 
                window), by = 1))
        } else {
            win <- queryHits(over)
        }
        
        subASMtuple <- as.data.frame(ASMtuple[win, ])
        submeth <- as.data.frame(meth[win, ])
        subASMtuple$pos <- submeth$pos <- tupgr$midpt[win]
        subASMtuple_long <- reshape2::melt(subASMtuple, id.vars = "pos", 
            measure.vars = colnames(ASMtuple))
        subASMtuple_long$score <- "ASMtuple"
        
        # Marginal meth
        submeth_long <- reshape2::melt(submeth, id.vars = "pos", 
            measure.vars = colnames(meth))
        submeth_long$score <- "meth"
        
    }
    
    if ((!is.null(derASM)) & (!is.null(ASM))) {
        message("Using both scores")
        full_long <- rbind(subASMtuple_long, submeth_long, subASMsnp_long, 
            subref_long, subalt_long)
        full_long$score <- factor(full_long$score, levels = c("ASMtuple", 
            "meth", "ASMsnp", "REF:meth", "ALT:meth"))
        full_long$group <- colData(derASM)$group[match(full_long$variable, 
            colData(derASM)$samples)]
    }
    
    if (is.null(ASM) & !is.null(derASM)) {
        message("Using ASMsnp score")
        full_long <- rbind(subASMsnp_long, subref_long, subalt_long)
        full_long$score <- factor(full_long$score, levels = c("ASMsnp", 
            "REF:meth", "ALT:meth"))
        full_long$group <- colData(derASM)$group[match(full_long$variable, 
            colData(derASM)$samples)]
    }
    
    if (is.null(derASM) & !is.null(ASM)) {
        message("Using ASMtuple score")
        full_long <- rbind(subASMtuple_long, submeth_long)
        full_long$group <- colData(ASM)$group[match(full_long$variable, 
            colData(ASM)$samples)]
    }
    
    # To draw rectangle on region
    forect <- data.frame(xmin = start(dame), xmax = end(dame), 
        ymin = 0, ymax = Inf)
    
    # plot scores
    m1 <- ggplot(data = full_long) + geom_line(mapping = aes_(x = ~pos, 
        y = ~value, group = ~variable, color = ~group), alpha = 0.5) + 
        geom_point(mapping = aes_(x = ~pos, y = ~value, group = ~variable, 
            color = ~group)) + geom_rect(data = forect, aes_(xmin = ~xmin, 
        xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), alpha = 0.1) + 
        facet_grid(score ~ ., scales = "free_y") + theme_bw() + 
        xlab("position") + ggtitle("Scores")
    
    cord <- ggplot_build(m1)$layout$panel_scales_x[[1]]$range$range
    
    # plot SNP table
    if (plotSNP) {
        m2 <- ggplot(data = subSNP_long) + 
            geom_point(aes_(x = ~pos, y = 1, group = ~variable, 
                        color = ~snp.pos), shape = 8, size = 1, 
                        fill = "white", stroke = 1) +
            facet_grid(variable ~ .) + 
            theme_bw() + 
            theme(panel.spacing = unit(0, "lines"), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()) + 
            coord_cartesian(xlim = cord) + 
            ylab("") + 
            ggtitle("SNPs - for ASMsnp") + 
            xlab("position")
        
        if (!is.null(colvec)) {
            m2 <- m2 + scale_color_manual(values = colvec)
            m1 <- m1 + scale_color_manual(values = colvec)
        }
        
        p <- cowplot::plot_grid(m1, m2, ncol = 1, nrow = 2, rel_heights = c(3, 
            1), align = "v")
    } else {
        
        if (!is.null(colvec)) {
            p <- m1 + scale_color_manual(values = colvec)
        } else {
            p <- m1
        }
    }
    
    return(p)
}
