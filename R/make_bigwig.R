#' Build bigwig files for visualization in IGV
#'
#' Builds a bigwig file for each sample, to visualize the ASM score and the derived ASM from
#' the SummarizedExperiment objects. Requires bedGraphToBigWig to be externally installed.
#'
#' @param sample Name of sample (as in colnames from SummarizedExperiment)
#' @param score.obj SummarizedExperiment with der.ASM or asm in assays(). Ideally filtered by coverage
#' previously
#' @param folder Destination folder for bigiwig file
#' @param chromsizes.file hg19.chrom.sizes.mod or any other chrom.sizes file as in (put link to get this file)
#'
#' @return A bigwig file for the sample chosen.
#'
#' @export
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

  #make bedgraph
  scorename <- names(SummarizedExperiment::assays(score.obj))[1]
  write.table(bg, file = sprintf("%s/%s.%s.bedgraph", folder, scorename, sample),
              sep = "\t", row.names = F, col.names = F, quote = F)

  #sort file
  cmd1 <- sprintf("sort -k1,1 -k2,2n %s/%s.%s.bedgraph > %s/%s.%s.sorted.bedgraph",
                  folder, scorename, sample, folder, scorename, sample)
  cat(cmd1, "\n")
  system(cmd1)

  #then make bigwig
  cmd2 <- sprintf("bedGraphToBigWig %s/%s.%s.sorted.bedgraph %s %s/%s.%s.bw",
                  folder, scorename, sample, chromsizes.file, folder, scorename, sample)
  cat(cmd2, "\n")
  system(cmd2)
}

#apply to all samples in SumExp
#sapply(colnames(ASM_score_matrix), make_bigwig, score.obj = ASM_score_matrix, folder = , chromsizes.file = )
