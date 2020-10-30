context("test-derivedASM_dames")
library(SummarizedExperiment)
library(limma)

data(extractbams_output)
grp <- factor(c(rep("CRC",4),rep("NORM",4)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)

derASM <- calc_derivedasm(extractbams_output, cores = 1, verbose = FALSE)
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) == 8 #filt to avoid warnings
derASM <- derASM[filt,]

test_that("end to end calc_derivedasm", {
  expect_s4_class(derASM, "RangedSummarizedExperiment")
  expect_length(colnames(derASM), 8)
})


test_that("all assays created", {
  expect_length(assays(derASM), 7)
})

# test_that("end to end tstat calc", {
#   derASMt <- get_tstats(derASM, mod, verbose = FALSE)
#   expect_s4_class(derASMt, "RangedSummarizedExperiment")
#   expect_s4_class(rowData(derASMt), "DataFrame")
# })

test_that("end to end find_dames", {
  dames <- find_dames(derASM, mod, verbose = FALSE)
  expect_is(dames, "data.frame")
})
