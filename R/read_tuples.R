#' Read in list of methtuple files
#'
#' This function reads in a list of files obtained from the methtuple tool. It
#' filters out tuples based on the set minimum coverage (min_cov) and the
#' maximum allowed distance (maxGap) between two genomic positions in a tuple.
#'
#'
#' @param files List of methtuple files.
#' @param sampleNames Names of files in the list.
#' @param minCoverage The minimum coverage per tuple. Tuples with a coverage <
#'   minCoverage are filtered out. Default = 2.
#' @param verbose If the function should be verbose.
#' @param maxGap The maximum allowed distance between two positions in a tuple.
#'   Only distances that are <= maxGap are kept. Default = 150 base pairs.
#'
#' @return A list of data frames, where each data frame corresponds to one file.
#' @examples
#' DATA_PATH_DIR <- system.file('extdata', '.', package = 'DAMEfinder')
#' get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
#'
#' tuple_files <- list.files(DATA_PATH_DIR, '.tsv.gz')
#' tuple_files <- get_data_path(tuple_files)
#' ASM <- read_tuples(tuple_files, c('CRC1', 'NORM1'))
#' 
#' @export
#'
#' 
read_tuples <- function(files, sampleNames, minCoverage = 2, 
    maxGap = 20, verbose = TRUE) {
    pos2 <- pos1 <- NULL
    MM <- UU <- MU <- UM <- NULL
    
    # read in methtuple files
    methtuple_cols <- readr::cols(chr = readr::col_character(), 
        strand = readr::col_character(), pos1 = readr::col_integer(), 
        pos2 = readr::col_integer(), MM = readr::col_integer(), 
        MU = readr::col_integer(), UM = readr::col_integer(), 
        UU = readr::col_integer())
    
    names(files) <- sampleNames
    
    methtuple_list <- lapply(files, function(u) {
        if (verbose) 
            message("Reading ", u)
        readr::read_tsv(u, col_types = methtuple_cols, progress = FALSE)
    })
    
    # filter out low coverage and maxGap
    if (verbose) 
        message("Filtering and sorting: ", appendLF = FALSE)
    methtuple_list <- lapply(methtuple_list, function(df) {
        if (verbose) 
            message(".", appendLF = FALSE)
        df$cov <- with(df, MM + UU + UM + MU)
        df$inter_dist <- with(df, pos2 - pos1)
        w <- df$cov >= minCoverage & df$inter_dist <= maxGap
        df <- df[w, ]
        df[order(df$chr, df$pos1, df$pos2), ]
    })
    if (verbose) 
        message(" done.")
    
    methtuple_list
}
