#' Build bigwig files for visualization in IGV
#'
#' Creates a \code{bigwig} file with ASM score for each sample. Requires \code{bedGraphToBigWig}. 
#' installed.
#' 
#' The \code{bedGraphToBigWig} script can be downloaded from 
#' [here](https://github.com/ENCODE-DCC/kentUtils/tree/master/bin/linux.x86_64). The \code{chrom.sizes}
#' file can be downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/), or a
#' similar path goldenPath. 
#'
#' @param sample Name of sample (as in colnames from \code{SummarizedExperiment})
#' @param score.obj \code{RangedSummarizedExperiment}, ideally filtered by coverage from 
#' \code{\link{calc_asm}} or \code{\link{calc_derivedasm}}.
#' @param folder Destination folder for bigwig file.
#' @param chromsizes.file hg19.chrom.sizes or any other file similar specifying the chromosome sizes.
#'
#' @return A bigwig file for the sample chosen.
#'
#' @export
#' @importFrom utils write.table
#'
#' @examples
make_bigwig <- function(sample, score.obj, folder, chromsizes.file){

  if(is.null(dim(SummarizedExperiment::assays(score.obj)[["asm"]]))){
    asm <- SummarizedExperiment::assays(score.obj)[["der.ASM"]]
    midpt <- BiocGenerics::start(score.obj)
  } else {
    asm <- SummarizedExperiment::assays(score.obj)[["asm"]]
    midpt <- S4Vectors::mcols(score.obj)$midpt
  }

  #make dataframe with midpoint and score for file
  bg <- data.frame(chr = GenomeInfoDb::seqnames(score.obj),
                   start = as.integer(midpt),
                   end = as.integer(midpt),
                   score = abs(asm[,grep(sample, colnames(asm))]),
                   stringsAsFactors = F,
                   row.names = NULL)

  #remove NAs (this is just to double check, the user should input a filtered matrix)
  bg <- bg[!is.na(asm[,grep(sample, colnames(asm))]),]

  #make bigwig
  #rtracklayer::export(object = GRasm, con = "bigwigs/CRC1.asmScore.bw", format = "bw") #bad function
  
  if(folder == "."){
    folder = ""
  } else {folder = paste0(folder,"/")}

  #make bedgraph
  scorename <- names(SummarizedExperiment::assays(score.obj))[1]
  write.table(bg, file = sprintf("%s%s.%s.bedgraph", folder, scorename, sample),
              sep = "\t", row.names = F, col.names = F, quote = F)

  #sort file
  cmd1 <- sprintf("sort -k1,1 -k2,2n %s%s.%s.bedgraph > %s%s.%s.sorted.bedgraph",
                  folder, scorename, sample, folder, scorename, sample)
  cat(cmd1, "\n")
  system(cmd1)

  #then make bigwig
  cmd2 <- sprintf("bedGraphToBigWig %s%s.%s.sorted.bedgraph %s %s%s.%s.bw",
                  folder, scorename, sample, chromsizes.file, folder, scorename, sample)
  cat(cmd2, "\n")
  system(cmd2)
}

#apply to all samples in SumExp
#sapply(colnames(ASM_score_matrix), make_bigwig, score.obj = ASM_score_matrix, folder = , chromsizes.file = )
