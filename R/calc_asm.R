#' Calculate ASM Score
#'
#' This function takes in a list of samples resulting from the read_tuples
#' function and returns a SummarizedExperiment of Allele-Specific Methylation
#' (ASM) scores, where each row is a tuple and each column is a sample.
#'
#' Calculates ASM score for a list of samples in the output format of the result
#' of read_tuples This functions uses the following other functions: process,
#' calcScore, calcWeight.
#'
#' @param sample_list List of samples returned from \code{\link{read_tuples}}
#' @param beta The beta parameter used to calculate the weight in the ASM score.
#'   \code{link{calc_weight}} uses this parameter for \code{\link{pbeta}}.
#'   Default = 0.5.
#' @param verbose If the function should be verbose. Default = TRUE.
#' @param a The distance from 0.5 allowed, where 0.5 is a perfect MM:UU balance
#'   for a tuple. In the default mode this value is set to 0.2, and we account
#'   for the instances where the balance is between 0.3 and 0.7.
#' @param transform Transform the calculated tuple ASM scores. We use the
#'   modulus square root function which outputs the square root, while
#'   preserving the original sign.
#'
#' @return A \code{SummarizedExperiment} of ASM scores where the rows are all
#'   the tuples and the columns the sample names.
#' @export
#'
#' @examples
calc_asm <- function(sample_list, beta=0.5, a=0.2, transform=modulus_sqrt, 
                     verbose=TRUE) {


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

  # get matrix of coverage, and tuple methylation
  if(verbose) message("Assembling coverage tables: ", appendLF = FALSE)
  coverage <- mapply( function(df,k){
    if(verbose) message(".", appendLF=FALSE)
    m <- match(key, k)
    df$cov[m]
  }, sample_list, all_keys)


  #MM
  MM <- mapply( function(df,k){
    if(verbose) message(".", appendLF=FALSE)
    m <- match(key, k)
    df$MM[m]
  }, sample_list, all_keys)

  #MU
  MU <- mapply( function(df,k){
    if(verbose) message(".", appendLF=FALSE)
    m <- match(key, k)
    df$MU[m]
  }, sample_list, all_keys)

  #UM
  UM <- mapply( function(df,k){
    if(verbose) message(".", appendLF=FALSE)
    m <- match(key, k)
    df$UM[m]
  }, sample_list, all_keys)

  #UU
  UU <- mapply( function(df,k){
    if(verbose) message(".", appendLF=FALSE)
    m <- match(key, k)
    df$UU[m]
  }, sample_list, all_keys)

  #Get ranges
  ss <- limma::strsplit2(rownames(asm),".",fixed=TRUE)
  gr <- GenomicRanges::GRanges(ss[,1],
                               IRanges::IRanges(as.numeric(ss[,2]),
                                                as.numeric(ss[,3])))
  names(gr) <- rownames(asm)
  gr$midpt <- floor((GenomicRanges::start(gr)+GenomicRanges::end(gr))/2)

  #Build object
  sa <- SummarizedExperiment::SummarizedExperiment(
    assays=S4Vectors::SimpleList(asm=abs(asm),
                                 cov = coverage,
                                 MM = MM,
                                 MU = MU,
                                 UM = UM,
                                 UU = UU),
    rowRanges=gr)

  if(verbose) message("Returning SummarizedExperiment with ",nrow(asm), 
                      " CpG pairs", appendLF = FALSE)
  o <- order(GenomeInfoDb::seqnames(sa),gr$midpt)
  sa[o]
}



#' Calculate score
#'
#' This function calculates the ASM score for every tuple in a given sample. The
#' ASM score is a multiplication of the log odds ratio by a weight that reflects
#' the extent of allele-specific methylation. This weight is obtained with the
#' \code{\link{calc_weight}} function.
#'
#' @param beta parameter for the \code{calc_weight} function. It's the alpha and
#'   beta values for the Beta function.
#' @param a parameter for the calc_weight function. The weight will be the
#'   probability that the MM/(MM+UU) ratio lies between 0.5-a and 0.5+a.
#' @param df data frame of a sample containing all information per tuple
#'   (MM,UU,UM and MU counts, as well as the log odds ratio per tuple) needed
#'   for the ASM score.
#' @details This function returns an allele-specific methylation (ASM) score for
#'   every given tuple in a sample. The ASM score is a product of the log odds
#'   ratio and a weight reflecting a measure of allele-specificity using the MM
#'   and UU counts.
#'
#'
#' @return The same object with an additional column for the ASM score.
#' @export
#'
#' @examples
#' # Get list of samples containing tuple count information
#' #sample_list <- read_tuples(files = file_list,
#' #sample_names = names(file_list))
#'
#' # Calculate log odds ratio per tuple for the first sample in the list
#' # sample1 <- calc_logodds(sample_list[[1]])
#'
#' # Calculate ASM score
#' # sample1 <- calc_score(sample1, beta=0.5, a=0.2)
#' 
calc_score <- function(df, beta = 0.5, a = 0.2) {
  weights <- calc_weight(df$MM, df$UU, beta=beta, a=a)
  df$asm_score <- df$logodds*weights #changed
  df
}

#' Calculate the log odds ratio
#'
#' This function calculates the log odds ratio for a CpG tuple:
#' \code{(MM*UU)/(UM*MU)}, where 'M' stands for methylated and 'U' for
#' unmethylated. 'MM' reflects the count for instances the CpG pair is
#' methylated at both positions. The higher the MM and UU counts for that CpG
#' pair, the higher the log odds ratio.
#'
#'
#' @param s A data frame that contains the MM,UU,UM, and MU counts for each CpG
#'   tuple for a particular sample. It is the resulting object of the
#'   \code{read_tuples}.
#' @param eps Count added to each of the MM,UU,UM and MU counts to avoid
#'   dividing by zero for example. The default is set to 1.
#'
#' @return The same object is returned with an additional column for the log
#'   odds ratio.
#' @export
#'
#' @examples
#' 
calc_logodds <- function(s, eps=1) {
  ratio <- with(s, ( (MM+eps)*(UU+eps) ) / ( (MU+eps)*(UM+eps) ) )
  s$logodds <- log10(ratio)
  s
}


#' Get Modulus Square Root
#'
#' Function to calculate signed square root (aka modulus square root).
#'
#' @param values Vector or matrix of ASM scores where each column is a sample.
#'   These values are transformed with a square root transformation that
#'   preserves the sign.
#'
#' @return Vector or matrix of transformed scores.
#' @export
#'
#' @examples
#' v <- c(-1,1.5,0,-2.1,4.3)
#' modulus_sqrt(v)
#'
#'
#' 
modulus_sqrt <- function(values) {
  sign(values)*sqrt(abs(values))
}


#' Calculate Weight for ASM Score
#'
#' This function calculates a weight which reflects MM to UU balance, where M
#' stands for methylated and U for unmethylated. Given the MM and UU counts for
#' a particular tuple, the weight is obtained using the \code{link{pbeta}}
#' function.
#'
#'
#' @param beta parameter for the beta distribution. In B(alpha,beta), we set
#'   alpha=beta=0.5 by default.
#' @param MM The read counts for where pos1 and pos2 of the tuple were both
#'   methylated.
#' @param UU The read counts for where pos1 and pos2 of the tuple were both
#'   unmethylated.
#' @param a parameter for how far from 0.5 we go as a measure of allele-specific
#'   methylation. The weight is the probability that the MM:(MM+UU) ratio is
#'   between 0.5-a and 0.5+a. The default is set to 0.2.
#'   
#' @details For a given tuple with MM and UU counts, the weight that reflects 
#' allele-scpecificity is calculated as follows:
#' 
#' - Prior:\deqn{p(\theta|\alpha,\beta) \sim Beta(\alpha,\beta),} where
#' \eqn{\theta = \frac{MM}{MM+UU}} and \eqn{\alpha = \beta = 0.5}. 
#' \eqn{p(\theta|\alpha,\beta)} represents our prior belief which is that 
#' tuples are either fully methylated or fully unmethylated, rather than 
#' allele-specifically methylated which is a much rarer event.
#' - Likelihood: \deqn{p(x|\alpha,\beta) \propto \theta^{MM}(1-\theta)^{UU},} 
#' where x is our observation (the MM and UU counts).
#' - Posterior:\deqn{p(\theta|x) \propto p(x|\theta)*p(\theta|\alpha,\beta)}
#' \deqn{p(\theta|x) \propto \theta^{MM-0.5}(1-\theta)^{UU-0.5},} 
#' where \eqn{\alpha = \beta = 0.5}. This posterior also follows a beta 
#' distribution \eqn{\sim Beta(\alpha'=MM+0.5, \beta'=UU+0.5)}
#' @md
#'
#' @return A number that reflects allele-specificity given MM and UU counts for
#'   a CpG pair. This is used as a weight that is multiplied by the log odds
#'   ratio to give the final ASM score of that tuple.
#'
#' @export
#'
#' @examples
#' weight1 <- calc_weight(MM=50, UU=50)
#' 0.9999716
#'
#' weight2 <- calc_weight(MM=20, UU=60)
#' 0.1646916
#'
#' weight3 <- calc_weight(MM=2, UU=3)
#' 0.6260661
#'
#' 
calc_weight <- function(MM, UU, beta=0.5, a=.2) {
  s1 <- beta+MM
  s2 <- beta+UU
  stats::pbeta(.5+a, shape1=s1, shape2=s2)-stats::pbeta(.5-a, 
                                                        shape1=s1, shape2=s2)
}
