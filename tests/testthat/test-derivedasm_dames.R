context("test-derivedASM_dames")
library(SummarizedExperiment)
library(limma)

data(splitbams_output)
grp <- factor(c(rep("CRC",4),rep("NORM",4)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)

test_that("end to end calc_derivedasm", {
  derASM <- calc_derivedasm(splitbams_output, cores = 1, verbose = FALSE)
  expect_s4_class(derASM, "RangedSummarizedExperiment")
  expect_length(colnames(derASM), 8)
})


derASM <- calc_derivedasm(splitbams_output, cores = 1, verbose = FALSE)
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) == 8 #filt to avoid warnings
derASM <- derASM[filt,]

test_that("all assays created", {
  expect_length(assays(derASM), 6)
})

test_that("end to end tstat calc", {
  derASMt <- get_tstats(derASM, mod, verbose = FALSE)
  expect_s4_class(derASMt, "RangedSummarizedExperiment")
  expect_s4_class(rowData(derASMt), "DataFrame")
  #expect_length(colnames(rowData(derASMt)), 3)
  #expect_equal(sum(is.na(rowData(derASMt)$smooth_tstat)), 0)
})

#TODO: test different lmfit methods and pvalAssign

test_that("end to end find_dames", {
  dames <- find_dames(derASM, mod, minNum = 2, minInSpan = 2, verbose = FALSE)
  expect_is(dames, "data.frame")
  #expect_equal(dim(dames)[1], 2)
})
