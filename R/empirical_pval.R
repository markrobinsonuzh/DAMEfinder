#' Calculate empirical region-level p-value
#'
#' This function permutes the coefficient of interest and re-runs
#' \code{\link{get_tstats}} and \code{\link{regionFinder}} for each permutation.
#' Code for permutations copied from the \code{dmrseq} function from the package
#' of the same name.
#'
#' @param presa SExperiment output from \code{calc_derivedasm} or
#'   \code{calc_asm}.
#' @param design design matrix.
#' @param rforiginal data.frame of DAMEs calculated with original design.
#' @param coeff Coefficient of interest to permute.
#' @param cont same as in \code{get_tstats}.
#' @param smooth Boolean.
#' @param maxPerms Maximum possible permutations generated. Default = 10.
#' @param Q Quantile for cuttof.
#' @param maxGap Same as other functions in the package.
#' @param method lmFit method.
#' @param ... Passed to \code{get_tstats} and then to \code{loessByCluster}.
#' @return Vector of empirical p-values.
#' @keywords internal
#' 
empirical_pval <- function(presa, design, rforiginal, coeff, 
    cont, smooth, maxPerms = 10, Q, maxGap, method, ...) {
    
    sampleSize <- table(design[, coeff])
    
    # Since we allow the user to choose the coeff from the
    # design.matrix, I will never have more than 2 coefs, as
    # dmrseq does
    
    # limit possible number of perms!
    if (length(unique(design[, coeff])) == 2 && choose(nrow(design), 
        min(sampleSize)) < 5e+05) {
        
        perms <- utils::combn(seq(1, nrow(design)), min(sampleSize))
        # Remove redundant permutations (if balanced)
        if (length(unique(table(design[, coeff]))) == 1) {
            perms <- perms[, seq_len(ncol(perms)/2)]
        }
        
        # restrict to unique permutations that don't include any
        # groups consisting of all identical conditions
        rmv <- NULL
        for (p in seq_len(ncol(perms))) {
            if (length(unique(design[perms[, p], coeff])) == 
                1) {
                rmv <- c(rmv, p)
            }
        }
        
        if (length(rmv) > 0) 
            perms <- perms[, -rmv]
        
        # subsample permutations based on similarity to original
        # partition gives preference to those with the least
        # similarity
        if (maxPerms < ncol(perms)) {
            similarity <- apply(perms, 2, function(x) {
                max(table(design[x, coeff]))
            })
            perms.all <- perms
            perms <- NULL
            levs <- sort(unique(similarity))
            l <- 1
            num <- 0
            while (!(num == maxPerms) && l <= length(levs)) {
                keep <- sample(which(similarity == levs[l]), 
                min(maxPerms - num, sum(similarity == levs[l])))
                perms <- cbind(perms, perms.all[, keep])
                l <- l + 1
                num <- ncol(perms)
            }
        }
    } else message("Too many samples!")
    
    # Detect permuted dames
    message("Generating ", ncol(perms), " permutations", appendLF = TRUE)
    
    areas <- apply(perms, 2, function(i) {
        
        # message('Permutation ', num_perms[i], appendLF = TRUE)
        reorder <- i
        designr <- design
        
        if (length(unique(design[, coeff])) == 2 && !nrow(perms) == 
            nrow(designr)) {
            designr[, coeff] <- 0
            designr[reorder, coeff] <- 1
        } else {
            designr[, coeff] <- designr[reorder, coeff]
        }
        
        sa_perm <- get_tstats(presa, design = designr, coef = coeff, 
            contrast = cont, maxGap = maxGap, smooth = smooth, 
            verbose = TRUE, method = method, ...)
        
        # choose smoothed if true
        if (smooth) {
            sm_tstat <- S4Vectors::mcols(sa_perm)$smooth_tstat
        } else {
            sm_tstat <- S4Vectors::mcols(sa_perm)$tstat
        }
        
        # choose position to find regions
        if (names(assays(sa_perm))[1] == "asm") {
            midpt <- S4Vectors::mcols(sa_perm)$midpt
        } else {
            midpt <- BiocGenerics::start(sa_perm)
        }
        
        K <- stats::quantile(abs(sm_tstat), Q, na.rm = TRUE)
        permrf <- bumphunter::regionFinder(x = sm_tstat, 
            chr = as.character(GenomeInfoDb::seqnames(sa_perm)), 
            pos = midpt, cluster = S4Vectors::mcols(sa_perm)$cluster, 
            cutoff = K, maxGap = maxGap, verbose = FALSE)
        return(abs(permrf$area))
    })
    
    # if( put a condition to check if areas actually has
    # something
    all_areas <- sort(unlist(areas))
    total_areas <- length(all_areas)
    
    pvalEmp <- vapply(rforiginal$area, function(a) {
        pperm <- (sum(all_areas > abs(a)) + 1)/(total_areas + 
            1)
        return(pperm)
    }, FUN.VALUE = double(1))
    
    return(pvalEmp)
}
