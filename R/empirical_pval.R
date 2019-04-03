#' Calculate emprical-region-level p-value
#'
#' This function performs permutations bla bla
#' Code for permutations copied from the \code{dmrseq} function
#' from the package of the same name.
#'
#' @param pre_sa
#' @param design 
#' @param rforiginal
#' @param coeff
#' @param smooth
#' @param maxPerms
#' @param K
#' @param maxGap
#' @param ... 
#' @return Vector of summarized pvals
#'
#' @md
#' @export
#' 
#' 
#' 
#' 
empirical_pval <- function(pre_sa, design, rforiginal, coeff, smooth, maxPerms = 10, K, maxGap, ...){
  
  sampleSize <- table(design[, coeff])
  
  if(length(unique(design[, coeff])) == 2 && 
      
      #limit possible number of perms!
      choose(nrow(design), min(sampleSize)) < 5e5 ){
  
    perms <- combn(seq(1, nrow(design)), min(sampleSize))
    
    # Remove redundant permutations (if balanced)
    if(length(unique(table(design[,coeff]))) == 1){
      perms <- perms[, seq_len(ncol(perms)/2)]
    }
    
    # restrict to unique permutations that don't include any 
    # groups consisting of all identical conditions
    rmv <- NULL
    for(p in seq_len(ncol(perms))){
      if(length(unique(design[perms[,p],coeff])) == 1){
        rmv <- c(rmv, p)
      }
    }
    
    if(length(rmv) > 0 ) perms <- perms[,-rmv]
    
    # subsample permutations based on similarity to original partition
    # gives preference to those with the least similarity
    if(maxPerms < ncol(perms)){
      similarity <- apply(perms, 2, function(x) {
        max(table(design[x,coeff]))
      })
      perms.all <- perms
      perms <- NULL
      levs <- sort(unique(similarity))
      l <- 1
      num <- 0
      while(!(num == maxPerms) && l <= length(levs)) {
        keep <- sample(which(similarity == levs[l]), 
                       min(maxPerms-num, sum(similarity == levs[l])) )
        perms <- cbind(perms, perms.all[,keep])
        l <- l + 1
        num <- ncol(perms) 
      }
    }
  } else {
    # Next consider a multilevel, or continuous covariate where the
    # covariate will be permuted in an unrestricted manner
 
       perms <- as.matrix(seq_len(nrow(design)))
    
    for (p in seq_len(maxPerms)) {
      tries <- 0
      candidate <- sample(seq_len(nrow(design)), nrow(design))
      # check that the permutation is not a duplicate, and not 
      # equal to the original
      while ((sum(apply(perms, 2, function(x) 
        all.equal(x, candidate)) == TRUE) > 0 || 
        sum(apply(perms, 2, function(x) 
          all.equal(x, rev(candidate))) == TRUE) > 0) &&
        tries <= 20) {
        candidate <- sample(seq(seq_len(nrow(design))), nrow(design))
        tries <- tries + 1
      }
      # save the permutation to the permutation matrix
      if (tries <= 20){
        perms <- cbind(perms, candidate)
      }
    }
    perms <- perms[,-1] # remove original
  }
  
  #Detect permuted dames
  message("Generating permutations")
  areas <- apply(perms, 2, function(i){
    
    reorder <- i
    designr <- design
    
    if(length(unique(design[, coeff])) == 2 && 
        !nrow(perms) == nrow(designr)){
      designr[,coeff] <- 0
      designr[reorder, coeff] <- 1
    } else {
      designr[, coeff] <- designr[reorder, coeff]
    }
    
    sa_perm <- get_tstats(pre_sa, 
                          design = designr,
                          coef = coeff,
                          maxGap = maxGap,
                          smooth  = smooth,
                          verbose = FALSE,
                          ...)

    # choose smoothed if true
    if(smooth){
      sm_tstat <- S4Vectors::mcols(sa_perm)$smooth_tstat
    } else {
      sm_tstat <- S4Vectors::mcols(sa_perm)$tstat
    }
    
    #choose position to find regions
    if(names(assays(sa_perm))[1] == "asm"){
      midpt <- S4Vectors::mcols(sa_perm)$midpt
    } else {
      midpt <- BiocGenerics::start(sa_perm)
    }

    rf <- bumphunter::regionFinder(x = sm_tstat,
                                   chr = as.character(GenomeInfoDb::seqnames(sa_perm)),
                                   pos = midpt,
                                   cluster = S4Vectors::mcols(sa_perm)$cluster,
                                   cutoff = K, #use same K always
                                   maxGap = maxGap,
                                   verbose = FALSE)

    rf$area
  })

  all_areas <- sort(unlist(areas))
  total_areas <- length(all_areas)

  rf$pvalEmp <- sapply(rf$area, function(a){
    pperm <- (sum(all_areas > a) + 1) / (total_areas + 1)
  })
  
  rf$FDR <- p.adjust(pval, method = "BH")
  return(rf)
}

### old code

#####   
#tot_length <- length(treat) + length(control)
#combs <- utils::combn(tot_length, length(control))


# #remove the true constrats
# rem <- apply(combs, 2, function(i){
#   all.equal(i,control)
# })
# rem <- which(rem == T)
# if(!identical(rem, integer(0))) combs <- combs[,-rem]
