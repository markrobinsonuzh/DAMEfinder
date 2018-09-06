#' Read in list of methtuple files
#'
#' This function reads in a list of files obtained from the methtuple tool. It filters out tuples
#' based on the set minimum coverage (min_cov) and the maximum allowed distance (max_gap) between
#' two genomic positions in a tuple.
#'
#'
#' @param files List of methtuple files
#' @param sample_names Names of files in the list.
#' @param min_coverage The minimum coverage per tuple. Tuples with a coverage < min_coverage are
#' filtered out. Default = 2.
#' @param verbose If the function should be verbose.
#' @param max_gap The maximum allowed distance between two positions in a tuple. Only distances that
#' are <= max_gap are kept. Default = 150 base pairs.
#'
#' @return A list of data frames, where each data frame corresponds to one file.
#' @export
#'
#' @examples
read_tuples <- function(files, sample_names, min_coverage=2, max_gap=20, verbose=TRUE ) {
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

  names(files) <- sample_names

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


#' Remove tuples overlapping with SNPs TODO: FIX this function to remove tuples not close to a SNP
#'
#' @param df Data frame for a sample resulting from the read_tuples function where each row is a tuple,
#' and the columns indicates the tuple positions on the genome, the counts for the different
#' methylation states, the coverage, and the distance between the two positions in a tuple.
#' @param snp_key Vector of SNPs where each element represents the SNP positions as a character in
#' the form of "chr.SNP_position". For example a SNP on position 231 in chr12 would be represented
#' as "chr12.231".
#' @param verbose If the function should be verbose.
#'
#' @return The same data frame as the input but without the tuples overlapping with SNPs.
#' @export
#'
#' @examples
remove_tuples <- function(df, snp_key=NULL, verbose=TRUE) {
  key_1 <- paste0(df$chr,".", df$pos1)
  df <- df[!(key_1 %in% snp_key),]

  key_2 <- paste0(df$chr,".", df$pos2)
  df <- df[!(key_2 %in% snp_key),]

  df
}
