context("test-tupleasm_dames")

library(SummarizedExperiment)


DATA_PATH_DIR <- "../../inst/extdata"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

tuple_files <- list.files(DATA_PATH_DIR, ".tsv.gz")
tuple_files <- get_data_path(tuple_files)

test_that("end to end read_tuples", {
  ASM <- read_tuples(tuple_files, c("CRC1", "CRC2", "CRC3", "NORM1", "NORM3"))
  expect_is(ASM, "list")
  expect_length(ASM, 5)
})

test_that("end to end calc_asm", {
  ASMscore <- calc_asm(ASM)
  expect_s4_class(ASMscore, "RangedSummarizedExperiment")
  expect_length(colnames(ASMscore), 5)
  expect_length(assays(ASMscore), 6)
})

test_that("end to end tstat calc", {
  ASMt <- get_tstats(ASMscore, 4:5, 1:3)
  expect_s4_class(ASMt, "RangedSummarizedExperiment")
})

test_that("end to end find_dames", {
  dames <- find_dames(ASMt, verbose = F)
  expect_is(dames, "data.frame")
  expect_equal(dim(dames)[1], 19)
})