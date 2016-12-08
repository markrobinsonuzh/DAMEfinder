## DAMEfinder functions

# input: list of methtuple files, names of samples (optional)
# output: list of data frames of samples containing tuple information remaining after filtering,
# ordered by chr, pos1, pos2


#' Read in list of methtuple files and filter out according to set min_cov and max_gap
#'
#' @param files
#' @param sample_names
#' @param min_coverage
#' @param verbose 
#' @param max_gap
#'
#' @return
#' @export
#'
#' @examples
read_tuples <- function(files, sample_names=seq.int(1:length(files)), min_coverage=10, max_gap=150, verbose=TRUE ) {
  # read in methtuple files
  
  methtuple_cols <- readr::cols(
    chr = readr::col_character(),
    strand = readr::col_character(),
    pos1 = readr::col_integer(),
    pos2 = readr::col_integer(),
    MM = readr::col_integer(),
    MU = readr::col_integer(),
    UM = readr::col_integer(),
    UU = readr::col_integer()
  )
    
  methtuple_list <- lapply(files, function(u) {
    if(verbose) message("Reading ",u)
    readr::read_tsv(u, col_types=methtuple_cols, progress = FALSE)
  })

  # filter out low coverage and maxGap
  if(verbose) message("Filtering and sorting: ", appendLF = FALSE)
  methtuple_list <- lapply(methtuple_list, function(df) {
    if(verbose) message(".",appendLF = FALSE)
    df$cov <- with(df, MM + UU + UM + MU)
    df$inter_dist <- with(df, pos2-pos1)
    w <- df$cov >= min_coverage & df$inter_dist <= max_gap
    df <- df[w,]
    df[order(df$chr,df$pos1,df$pos2),]
  })
  if(verbose) message(" done.")

  methtuple_list
}


#' Title
#'
#' @param df 
#' @param snp_key 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
remove_snps <- function(df, snp_key=NULL, verbose=TRUE) {
  key_1 <- paste0(df$chr,".", df$pos1)
  df <- df[!(key_1 %in% snp_key),]
  
  key_2 <- paste0(df$chr,".", df$pos2)
  df <- df[!(key_2 %in% snp_key),]

  df
}