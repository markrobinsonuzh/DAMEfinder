
# B # Calculate ASM score for a list of samples in the output format of the result of read_tuples
# This functions uses the following other functions: process, calcScore, calcWeight

#' Title
#'
#' @param sample_list
#' @param beta
#' @param a
#' @param return_matrix
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calc_asm <- function(sample_list, beta=0.5, a=0.2, return_matrix=TRUE, 
                     order_by_midpoint=TRUE, verbose=TRUE) {

  if(verbose) message("Calculating log odds.")
  sample_list <- lapply(sample_list, calc_logodds)

  if(verbose) message("Calculating ASM score: ", appendLF=FALSE)
  sample_list <- lapply(sample_list, function(u) {
    if(verbose) message(".", appendLF=FALSE)
    calc_score(u, beta=beta, a=a)
  })
  if(verbose) message(" done.")
  
  if (return_matrix) {

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
    
    if( order_by_midpoint) {
      asm <- order_pos_by_median(asm)
    }
    
    return(asm)
  }
  sample_list
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
calc_score <- function(df, beta, a) {
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



#' Transform the ASM scores
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
transform_scores <- function(x) {
  modulus_sqrt(x)
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
  pbeta(.5+a, shape1=s1, shape2=s2)-pbeta(.5-a, shape1=s1, shape2=s2)
}



#' Title
#'
#' @param score_matrix
#'
#' @return
#' @export
#'
#' @examples
order_pos_by_median <- function(score_matrix) {
  
  # set median of each tuple as the genomic position
  pos <- rownames(score_matrix)
  ss <- limma::strsplit2(pos, ".", fixed=TRUE)
  chr <- ss[,1]
  pos1 <- as.numeric(ss[,2])
  pos2 <- as.numeric(ss[,3])
  midpt <- (pos1+pos2)/2
  
  # sort the score matrix by median position (important for regionFinder and bumphunting)
  o <- order(chr, midpt)
  score_matrix[o,]
  
}
