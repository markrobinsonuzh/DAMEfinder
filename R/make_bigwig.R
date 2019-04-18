#' Build bigwig files for visualization in IGV
#'
#' Creates a \code{bigwig} file with ASM score for each sample. Requires
#' \code{bedGraphToBigWig}. installed.
#'
#' The \code{bedGraphToBigWig} script can be downloaded from
#' \url{https://github.com/ENCODE-DCC/kentUtils/tree/master/bin/linux.x86_64}.
#' The \code{chrom.sizes} file can be downloaded from
#' \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/}, or a similar
#' goldenPath.
#'
#' @param sample Name of sample (as in colnames from
#'   \code{SummarizedExperiment})
#' @param scoreObj \code{RangedSummarizedExperiment}, ideally filtered by
#'   coverage from \code{\link{calc_asm}} or \code{\link{calc_derivedasm}}.
#' @param folder Destination folder for bigwig file.
#' @param chromsizesFile hg19.chrom.sizes or any other file similar specifying
#'   the chromosome sizes.
#'
#' @return A bigwig file for the sample chosen.
#'
#' @importFrom utils write.table
#' @examples
#' #To apply to all samples in SummarizedExperiment
#' #sapply(colnames(ASM_score_matrix), make_bigwig, 
#' #scoreObj = ASM_score_matrix, folder = , chromsizesFile = )
#' @export
#'
make_bigwig <- function(sample, scoreObj, folder, chromsizesFile){

  if(is.null(dim(SummarizedExperiment::assays(scoreObj)[["asm"]]))){
    asm <- SummarizedExperiment::assays(scoreObj)[["der.ASM"]]
    midpt <- BiocGenerics::start(scoreObj)
  } else {
    asm <- SummarizedExperiment::assays(scoreObj)[["asm"]]
    midpt <- S4Vectors::mcols(scoreObj)$midpt
  }

  #make dataframe with midpoint and score for file
  bg <- data.frame(chr = GenomeInfoDb::seqnames(scoreObj),
                   start = as.integer(midpt),
                   end = as.integer(midpt),
                   score = abs(asm[,grep(sample, colnames(asm))]),
                   stringsAsFactors = F,
                   row.names = NULL)

  #remove NAs (this is just to double check, the user should input a filtered
  #matrix)
  bg <- bg[!is.na(asm[,grep(sample, colnames(asm))]),]
  
  if(folder == "."){
    folder = ""
  } else {folder = paste0(folder,"/")}

  #make bedgraph
  scorename <- names(SummarizedExperiment::assays(scoreObj))[1]
  write.table(bg, file = sprintf("%s%s.%s.bedgraph", folder, scorename, sample),
              sep = "\t", row.names = F, col.names = F, quote = F)

  #sort file
  cmd1 <- sprintf(
    "sort -k1,1 -k2,2n %s%s.%s.bedgraph > %s%s.%s.sorted.bedgraph",
    folder, scorename, sample, folder, scorename, sample)
  cat(cmd1, "\n")
  system(cmd1)

  #then make bigwig
  cmd2 <- sprintf("bedGraphToBigWig %s%s.%s.sorted.bedgraph %s %s%s.%s.bw",
                  folder, scorename, sample, chromsizesFile, folder, 
                  scorename, sample)
  cat(cmd2, "\n")
  system(cmd2)
}

