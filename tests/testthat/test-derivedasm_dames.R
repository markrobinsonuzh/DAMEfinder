context("test-derivedASM_dames")
library(DAMEfinder)
library(SummarizedExperiment)
library(limma)


DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

rds_files <- list.files(DATA_PATH_DIR, ".rds$")
rds_files <- get_data_path(rds_files)

test_that("end to end calc_derivedasm", {
  derASM <- calc_derivedasm(rds_files, cores = 1, verbose = T)
  expect_s4_class(derASM, "RangedSummarizedExperiment")
  expect_length(colnames(derASM), 8)
})

derASM <- calc_derivedasm(rds_files, cores = 1, verbose = T)

test_that("all assays created", {
  expect_length(assays(derASM), 6)
  expect_equal(sum(assays(derASM)[["der.ASM"]] > 1), 0)
})

test_that("end to end tstat calc", {
  derASMt <- get_tstats(derASM, 5:8, 1:4)
  expect_s4_class(derASMt, "RangedSummarizedExperiment")
})

#TODO: test different lmfit methods, and coefs

derASMt <- get_tstats(derASM, 5:8, 1:4, minNum = 2, minInSpan = 2)

test_that("all rowData fields created and matching", {
  expect_s4_class(rowData(derASMt), "DataFrame")
  expect_length(colnames(rowData(derASMt)), 3)
  expect_equal(sum(is.na(rowData(derASMt)$smooth_tstat)), 0)
})


test_that("end to end find_dames", {
  dames <- find_dames(derASMt, verbose = F)
  expect_is(dames, "data.frame")
  expect_equal(dim(dames)[1], 2)
})
