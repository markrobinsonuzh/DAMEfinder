context("test-dame_track")
library(GenomicRanges)

data(extractbams_output)

derASM <- calc_derivedasm(extractbams_output, cores = 1, verbose = TRUE)
DAME <- GRanges(19, IRanges(306443,310272))

SummarizedExperiment::colData(derASM)$group <- c(rep("CRC",4), rep("NORM", 4))

test_that("SNP dame track", {
  p <- dame_track(dame = DAME, derASM = derASM)
  expect_true(is.ggplot(p))
})

data(readtuples_output)
ASM <- calc_asm(readtuples_output)
SummarizedExperiment::colData(ASM)$group <- c(rep("CRC",3), rep("NORM", 2))
SummarizedExperiment::colData(ASM)$samples <- colnames(ASM)

test_that("dame track with both", {
  p <- dame_track(dame = DAME, derASM = derASM, ASM = ASM)
  expect_true(is.ggplot(p))
})

test_that("other options", {
  p <- dame_track(dame = DAME, derASM = derASM, ASM = ASM,
                  plotSNP = TRUE)
  expect_true(is.ggplot(p))
})
