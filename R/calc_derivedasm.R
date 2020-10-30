#' Calculate SNP-based ASM
#'
#' Combines all the \code{GRangeslist} generated in \code{\link{extract_bams}}
#' into a \code{\link{RangedSummarizedExperiment}} object, and calculates
#' SNP-based allele-specific methylation. 
#'
#' @param sampleList List of samples returned from \code{\link{extract_bams}}.
#' @param cores Number of cores to thread.
#' @param verbose If the function should be verbose.
#' @return \code{RangedSummarizedExperiment} containing in assays:
#'
#'   - der.ASM: matrix with SNP-based ASM 
#'   - snp.table: Matrix with SNP associated to the CpG site. 
#'   - ref.cov: Coverage of the 'reference' allele. 
#'   - alt.cov: Coevarage of the 'alternative' allele. 
#'   - ref.meth: Methylated reads from the 'reference' allele. 
#'   - alt.meth: Methylated reads from the 'alternative' allele.
#' @md
#'
#' @examples
#' data(extractbams_output)
#' derASM <- calc_derivedasm(extractbams_output[c(1,2)], cores = 1, 
#'    verbose = FALSE)
#' 
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors mcols<-
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @export

calc_derivedasm <- function(sampleList, cores = 1, verbose = TRUE) {
    
    # filter duplicated sites and choose the one with the highest
    # meth.diff
    allGR <- parallel::mclapply(sampleList, function(snp.table) {
        if (verbose) 
            message(".", appendLF = FALSE)
        
        # remove NULL values from list
        snp.table <- plyr::compact(snp.table)
        
        # Create a (actual) GRangesList, to be able to unlist
        GRlist <- GenomicRanges::GRangesList(snp.table)
        w <- unlist(GRlist)
        
        # calculate meth.diff (aka derASM)
        prop.alt <- mcols(w)$meth.alt/mcols(w)$cov.alt
        prop.ref <- mcols(w)$meth.ref/mcols(w)$cov.ref
        mcols(w)$meth.diff <- abs(prop.alt - prop.ref)
        
        # transform to keys
        full.keys <- paste0(seqnames(w), ".", start(w))
        
        # get unique ranges to loop
        unique.keys <- unique(sort(full.keys))
        
        # For duplicate sites, choose the one with the highest
        # meth.diff
        cols <- vapply(unique.keys, function(x) {
            # cols <- vapply(unGR, function(x){
            sHits <- as.integer(match(x, full.keys))
            
            if (length(sHits) <= 1) {
                return(methods::as(mcols(w)[sHits, ], "matrix"))
            } else {
                methdiffs <- mcols(w)$meth.diff[sHits]
                high <- as.integer(which.max(methdiffs))
                return(methods::as(mcols(w)[sHits[high], ], "matrix"))
            }
        }, character(6))
        
        # the vapply returns a tranversed matrix with everything as
        # characters, so I have to transform this.
        mcol <- data.frame(as.numeric(cols[1, ]), 
                           as.numeric(cols[2, ]), 
                           as.numeric(cols[3, ]), 
                           as.numeric(cols[4, ]), 
                           cols[5, ], 
                           as.numeric(cols[6, ]), stringsAsFactors = FALSE)
        
        colnames(mcol) <- colnames(mcols(w))
        ss <- limma::strsplit2(unique.keys, ".", fixed = TRUE)
        
        # return ordered GRanges
        unGR <- GRanges(ss[, 1], IRanges(as.integer(ss[, 2]), 
            width = 1))
        mcols(unGR) <- mcol
        
        # human ordering
        unGR <- GenomeInfoDb::sortSeqlevels(unGR)
        unGR <- sort(unGR)
        
        # The unGR has the filtered sites for this sample
        return(unGR)
    }, mc.cores = cores)
    
    # Try to extract from all the GRanges the unique sites by
    # generating a key
    all_keys <- lapply(allGR, function(u) {
        paste0(seqnames(u), ".", start(u))
    })
    
    key <- unique(unlist(all_keys))
    ss <- limma::strsplit2(key, ".", fixed = TRUE)
    keyGR <- GRanges(ss[, 1], IRanges(as.numeric(ss[, 2]), width = 1))
    
    # get matrix of scores across all samples
    if (verbose) 
        message("Summarizing scores")
    trueASM_table <- mapply(function(df, k) {
        m <- match(key, k)
        mcols(df)$meth.diff[m]
    }, allGR, all_keys)
    
    rownames(trueASM_table) <- key
    
    # Get matrices of coverages
    if (verbose) 
        message("Summarizing coverages")
    ref.cov_table <- mapply(function(df, k) {
        m <- match(key, k)
        return(mcols(df)$cov.ref[m])
    }, allGR, all_keys)
    
    alt.cov_table <- mapply(function(df, k) {
        m <- match(key, k)
        return(mcols(df)$cov.alt[m])
    }, allGR, all_keys)
    
    ref.meth_table <- mapply(function(df, k) {
        m <- match(key, k)
        return(mcols(df)$meth.ref[m])
    }, allGR, all_keys)
    
    alt.meth_table <- mapply(function(df, k) {
        m <- match(key, k)
        return(mcols(df)$meth.alt[m])
    }, allGR, all_keys)
    
    #Calculate stat ASM (prop test style)
    prop <- (ref.meth_table + alt.meth_table) / 
      (ref.cov_table + alt.cov_table)
    
    sigma <- sqrt(prop * (1 - prop) * ((1/ref.cov_table) + 
                                         (1/alt.cov_table)))
    
    zstat <- ((ref.meth_table/ref.cov_table) - 
                        (alt.meth_table/alt.cov_table)) / sigma
    
    # Get matrix for snp ID
    if (verbose) 
        message("Summarizing SNP info")
    snp_match_table <- mapply(function(df, k) {
        m <- match(key, k)
        mcols(df)$snp[m]
    }, allGR, all_keys)
    
    # length(snp_match_table)
    if (dim(snp_match_table)[1] != dim(trueASM_table)[1]) {
        stop("Tables contain different sizes")
    }
    
    if (dim(snp_match_table)[1] != dim(trueASM_table)[1]) {
        stop("Tables contain different sizes")
    } else {
        if (verbose) 
            message(sprintf("Returning %i CpG sites for %i samples", 
                dim(trueASM_table)[1], dim(trueASM_table)[2]))
    }
    
    ## Put all matrices into S4 vector
    derived_ASM_matrix <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(
        der.ASM = trueASM_table, 
        z.ASM = abs(zstat),                               
        snp.table = snp_match_table, 
        ref.cov = ref.cov_table, 
        alt.cov = alt.cov_table, 
        ref.meth = ref.meth_table, 
        alt.meth = alt.meth_table), 
        rowRanges = keyGR, 
        colData = S4Vectors::DataFrame(samples = names(allGR)))
    
    # last human sorting
    derived_ASM_matrix <- GenomeInfoDb::sortSeqlevels(derived_ASM_matrix)
    derived_ASM_matrix <- sort(derived_ASM_matrix)
    
    return(derived_ASM_matrix)
}

