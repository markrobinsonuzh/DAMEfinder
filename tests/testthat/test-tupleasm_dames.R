context("test-tupleasm_dames")

library(SummarizedExperiment)


DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

tuple_files <- list.files(DATA_PATH_DIR, ".tsv.gz")
tuple_files <- get_data_path(tuple_files)

data(readtuples_output)
grp <- factor(c(rep("CRC",3),rep("NORM",2)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)

test_that("end to end read_tuples", {
  ASM <- read_tuples(tuple_files, c("CRC1", "NORM1"))
  expect_is(ASM, "list")
  expect_length(ASM, 2)
})


ASMscore <- calc_asm(readtuples_output)
filt <- rowSums(!is.na(assay(ASMscore, "asm"))) == 5 #filt to avoid warnings
ASMscore <- ASMscore[filt,]

test_that("end to end calc_asm", {
  #ASMscore <- calc_asm(readtuples_output)
  expect_s4_class(ASMscore, "RangedSummarizedExperiment")
  expect_length(colnames(ASMscore), 5)
  expect_length(assays(ASMscore), 6)
})

# test_that("end to end tstat calc", {
#   ASMt <- get_tstats(ASMscore, mod, maxGap = 300, verbose = FALSE)
#   expect_s4_class(ASMt, "RangedSummarizedExperiment")
# })


test_that("end to end find_dames", {
  dames <- find_dames(ASMscore, mod, verbose = FALSE)
  expect_is(dames, "data.frame")
})

test_that("use permutations", {
  dames <- find_dames(ASMscore, mod, pvalAssign = "empirical", verbose = FALSE)
  expect_is(dames, "data.frame")
})
