context("test-derivedASM_dames")
library(SummarizedExperiment)
library(limma)


DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

sample_list <- try(readRDS(get_data_path("splitbams_output.rds")))

test_that("end to end calc_derivedasm", {
  derASM <- calc_derivedasm(sample_list, cores = 1, verbose = T)
  expect_s4_class(derASM, "RangedSummarizedExperiment")
  expect_length(colnames(derASM), 8)
})


derASM <- calc_derivedasm(sample_list, cores = 1, verbose = T)

test_that("all assays created", {
  expect_length(assays(derASM), 6)
})

test_that("end to end tstat calc", {
  filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) == 8 #filt to avoid warnings
  derASM <- derASM[filt,]
  derASMt <- get_tstats(derASM, 5:8, 1:4, verbose = F)
  expect_s4_class(derASMt, "RangedSummarizedExperiment")
  expect_s4_class(rowData(derASMt), "DataFrame")
  expect_length(colnames(rowData(derASMt)), 3)
  expect_equal(sum(is.na(rowData(derASMt)$smooth_tstat)), 0)
})

#TODO: test different lmfit methods, and coefs

test_that("end to end find_dames", {
  filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) == 8 #filt to avoid warnings
  derASM <- derASM[filt,]
  dames <- find_dames(derASM, 5:8, 1:4, minNum = 2, minInSpan = 2, verbose = F)
  expect_is(dames, "data.frame")
  expect_equal(dim(dames)[1], 2)
})
