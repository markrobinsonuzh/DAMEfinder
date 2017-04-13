
# B # Calculate ASM score for a list of samples in the output format of the result of read_tuples
# This functions uses the following other functions: process, calcScore, calcWeight

#' Title
#'
#' @param sample_list
#' @param beta
#' @param verbose 
#' @param a
#' @param transform 
#'
#' @return
#' @export
#'
#' @examples
calc_asm <- function(sample_list, beta=0.5, a=0.2, transform=modulus_sqrt, verbose=TRUE) {

  if(verbose) message("Calculating log odds.")
  sample_list <- lapply(sample_list, calc_logodds)

  if(verbose) message("Calculating ASM score: ", appendLF=FALSE)
  sample_list <- lapply(sample_list, function(u) {
    if(verbose) message(".", appendLF=FALSE)
    calc_score(u, beta=beta, a=a)
  })
  if(verbose) message(" done.")
  
  if(verbose) message("Creating position pair keys: ", appendLF = FALSE)
  # get key of unique tuples
  all_keys <- lapply(sample_list, function(u) {
    if(verbose) message(".", appendLF=FALSE)
    paste0(u$chr,'.',u$pos1, '.', u$pos2)
  })
  key <- unique(unlist(all_keys))
  if(verbose) message(" done.")
    
  # get matrix of ASM scores across all samples
  if(verbose) message("Assembling table: ", appendLF = FALSE)
  asm <- mapply( function(df,k){
    if(verbose) message(".", appendLF=FALSE)
    m <- match(key, k)
    df$asm_score[m]
  }, sample_list, all_keys)

  rownames(asm) <- key
  colnames(asm) <- names(sample_list)
  if(verbose) message(" done.")
    
  if(verbose) message("Transforming.")
  asm <- transform(asm)
  
  ss <- limma::strsplit2(rownames(asm),".",fixed=TRUE)
  gr <- GenomicRanges::GRanges(ss[,1], 
                               IRanges::IRanges(as.numeric(ss[,2]),
                                                as.numeric(ss[,3])))
  names(gr) <- rownames(asm)
  gr$midpt <- (GenomicRanges::start(gr)+GenomicRanges::end(gr))/2

  sa <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(asm=asm),
                                                   rowRanges=gr)  
  if(verbose) message("Returning SummarizedExperiment with ",nrow(asm), " CpG pairs", appendLF = FALSE)
  o <- order(GenomeInfoDb::seqnames(sa),gr$midpt)
  sa[o]
}



#' Get a vector of the posterior scores for all positions given beta and a
#'
#' @param beta
#' @param a
#' @param df
#'
#' @return
#' @export
#'
#' @examples
calc_score <- function(df, beta = 0.5, a = 0.2) {
  weights <- calc_weight(df$MM, df$UU, beta=beta, a=a)
  df$asm_score <- df$logodds*weights
  df
}

#' Calculate the odds ratio
#'
#' @param s
#'
#' @return
#' @export
#'
#' @examples
calc_logodds <- function(s, eps=1) {
  ratio <- with(s, ( (MM+eps)*(UU+eps) ) / ( (MU+eps)*(UM+eps) ) )
  s$logodds <- log10(ratio)
  s
}


#' Function to calculate signed square root (aka modulus square root)
#'
#' @param values
#'
#' @return
#' @export
#'
#' @examples
modulus_sqrt <- function(values) {
  sign(values)*sqrt(abs(values))
}


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
calc_weight <- function(MM, UU, beta=0.5, a=.2) {
  s1 <- beta+MM
  s2 <- beta+UU
  stats::pbeta(.5+a, shape1=s1, shape2=s2)-stats::pbeta(.5-a, shape1=s1, shape2=s2)
}
