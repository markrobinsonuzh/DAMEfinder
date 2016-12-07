## DAMEfinder functions


# A # Read in list of methtuple files and filter out according to set min_cov and max_gap.

# input: list of methtuple files, names of samples (optional)
# output: list of data frames of samples containing tuple information remaining after filtering,
# ordered by chr, pos1, pos2


#' Title
#'
#' @param files
#' @param sample_names
#' @param min_coverage
#' @param max_gap
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
read_tuples <- function(files, sample_names=seq.int(1:length(files)), min_coverage=10, max_gap=150, ... ) {
  # read in methtuple files
  methtuple_list <- parallel::mclapply(files, function(f){read.table(f, header=T)},
                                       mc.cores=min(length(files), ceiling(parallel::detectCores()/3)))

  # order each data frame by chr and pos1, add coverage column, and add tuple distance column
  methtuple_list_o <- lapply(methtuple_list, function(df) {
    df$cov <- with(df, MM + UU + UM + MU)
    df$inter_dist <- with(df, pos2-pos1)
    df[order(df[,"chr"], df[,"pos1"], df[,"pos2"]),]
  })

  # filter out low coverage and maxGap
  methtuple_list_o_f <- lapply(methtuple_list_o, function(df) {
    w <- which(df$cov >= min_coverage & df$inter_dist <= max_gap)
    return(df[w,])
  })

  return(methtuple_list_o_f)

}
