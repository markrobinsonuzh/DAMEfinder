context("test-dame_track_median")
library(GenomicRanges)


DAME <- GRanges(19, IRanges(306443,310272)) 
data(readtuples_output)
ASM <- calc_asm(readtuples_output)
SummarizedExperiment::colData(ASM)$group <- c(rep("CRC",3), rep("NORM", 2))
SummarizedExperiment::colData(ASM)$samples <- colnames(ASM)

test_that("dame track with ASM", {
  p <- dame_track_mean(dame = DAME, ASM = ASM)#derASM = derASM, )
  expect_true(is.ggplot(p))
})
