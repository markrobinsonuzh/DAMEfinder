#' Plot means per group of score tracks.
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
#' dame_track_mean(dame = DAME, 
#'                 ASM = ASM)
#' 
#' @export

dame_track_mean <- function(dame, window = 0, positions = 0, 
    derASM = NULL, ASM = NULL, colvec = NULL) {
    
    res_dame <- dame
    start(res_dame) <- start(dame) - positions
    end(res_dame) <- end(dame) + positions
    
    summarize_dat <- function(cond, SumExp, scorename, subtab, 
        pos) {
        idx <- colData(SumExp)$group == cond
        medians <- BiocGenerics::rowMeans(as.matrix(subtab)[, 
            idx], na.rm = TRUE)
        mins <- medians - (apply((as.matrix(subtab)[, idx]), 
            1, FUN = sd_rem))
        maxs <- medians + (apply((as.matrix(subtab)[, idx]), 
            1, FUN = sd_rem))
        return(data.frame(Means = medians, lower = mins, upper = maxs, 
            Group = cond, Score = scorename, pos = pos))
    }
    
    # custom se function
    sd_rem <- function(v) {
        if (sum(is.na(v)) == length(v)) {
            return(NA)
            stop()
        }
        v <- stats::na.omit(v)
        return(sqrt(stats::var(v)/length(v)))
    }
    
    if (!is.null(derASM)) {
        
        if (!all.equal(colData(derASM)$samples, colnames(derASM))) {
            stop("Sample names in colData() and colnames are different")
        }
        
        ASMsnp <- assay(derASM, "der.ASM")
        ref <- assay(derASM, "ref.meth")/assay(derASM, "ref.cov")
        alt <- assay(derASM, "alt.meth")/assay(derASM, "alt.cov")
        
        snpgr <- SummarizedExperiment::rowRanges(derASM)
        
        over <- GenomicRanges::findOverlaps(snpgr, res_dame)
        if (window != 0) {
            win <- c(seq(from = (queryHits(over)[1] - window), 
                to = (queryHits(over)[1] - 1), by = 1), queryHits(over), 
                seq(from = (utils::tail(queryHits(over), n = 1) + 1), 
                to = (utils::tail(queryHits(over), n = 1) + window), by = 1))
        } else {
            win <- queryHits(over)
        }
        
        # ASMsnp
        subASMsnp <- as.data.frame(ASMsnp[win, ])
        subref <- as.data.frame(ref[win, ])
        subalt <- as.data.frame(alt[win, ])
        substart <- start(snpgr)[win]
        
        ASMsnp_meds <- lapply(unique(colData(derASM)$group), 
            summarize_dat, SumExp = derASM, scorename = "ASMsnp", 
            subtab = subASMsnp, pos = substart)
        
        ref_meds <- lapply(unique(colData(derASM)$group), summarize_dat, 
            SumExp = derASM, scorename = "REF:meth", subtab = subref, 
            pos = substart)
        
        alt_meds <- lapply(unique(colData(derASM)$group), summarize_dat, 
            SumExp = derASM, scorename = "ALT:meth", subtab = subalt, 
            pos = substart)
        
        ASMsnp_meds_full <- do.call(rbind, ASMsnp_meds)
        ref_meds_full <- do.call(rbind, ref_meds)
        alt_meds_full <- do.call(rbind, alt_meds)
        
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
                seq(from = (utils::tail(queryHits(over), n = 1) + 1), 
                to = (utils::tail(queryHits(over), n = 1) + window), by = 1))
        } else {
            win <- queryHits(over)
        }
        
        subASMtuple <- as.data.frame(ASMtuple[win, ])
        submeth <- as.data.frame(meth[win, ])
        subtupger <- tupgr$midpt[win]
        
        ASMtuple_meds <- lapply(unique(colData(ASM)$group), summarize_dat, 
            SumExp = ASM, scorename = "ASMtuple", subtab = subASMtuple, 
            pos = subtupger)
        
        ASMtuple_meds_full <- do.call(rbind, ASMtuple_meds)
        
        meth_meds <- lapply(unique(colData(ASM)$group), summarize_dat, 
            SumExp = ASM, scorename = "meth", subtab = submeth, 
            pos = subtupger)
        
        meth_meds_full <- do.call(rbind, meth_meds)
    }
    
    
    if ((!is.null(derASM)) & (!is.null(ASM))) {
        message("Using both scores")
        full_long <- rbind(ASMtuple_meds_full, meth_meds_full, 
            ASMsnp_meds_full, ref_meds_full, alt_meds_full)
        full_long$Score <- factor(full_long$Score, levels = c("ASMtuple", 
            "meth", "ASMsnp", "REF:meth", "ALT:meth"))
    }
    
    if (is.null(ASM) & !is.null(derASM)) {
        message("Using ASMsnp score")
        full_long <- rbind(ASMsnp_meds_full, ref_meds_full, alt_meds_full)
        full_long$Score <- factor(full_long$Score, levels = c("ASMsnp", 
            "REF:meth", "ALT:meth"))
    }
    
    if (is.null(derASM) & !is.null(ASM)) {
        message("Using ASMtuple score")
        full_long <- rbind(ASMtuple_meds_full, meth_meds_full)
        
    }
    
    # To draw rectangle on region
    forect <- data.frame(xmin = start(dame), xmax = end(dame), 
        ymin = 0, ymax = Inf)
    
    # plot scores
    m1 <- ggplot(data = full_long) + geom_line(mapping = aes_(x = ~pos, 
        y = ~Means, color = ~Group)) + geom_ribbon(mapping = aes_(x = ~pos, 
        ymin = ~lower, ymax = ~upper, fill = ~Group), alpha = 0.2) + 
        geom_rect(data = forect, aes_(xmin = ~xmin, xmax = ~xmax, 
            ymin = ~ymin, ymax = ~ymax), alpha = 0.1) + facet_grid(Score ~ 
        ., scales = "free_y") + theme_bw() + xlab("position") + 
        ggtitle("Summary scores")
    
    if (!is.null(colvec)) {
        p <- m1 + scale_color_manual(values = colvec) + 
            scale_fill_manual(values = colvec)
    } else {
        p <- m1
    }
    
    return(p)
}
