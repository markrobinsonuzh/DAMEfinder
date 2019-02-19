context("test-split_reads")
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")
get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
bam.file <- get_data_path("NORM1_chr19_trim.bam")

snp <- GRanges(19, IRanges(388065, width = 1))

params <- ScanBamParam(
  tag = c("MD","XM","XR","XG"),
  which = snp)

alns.pairs <- readGAlignmentPairs(bam.file,
                                  use.names = TRUE,
                                  param = params)

alns <- unlist(alns.pairs)

test_that("end to end split_reads", {
  split <- splitReads(alns, "A", snp)
  expect_is(split, "list")
  expect_length(split$ref.reads, 72)
  expect_length(split$alt.reads, 14)
})

test_that("end to end get_MD", {
  a <- "2C0C3C7C2C7C1C0C1C1C1C2C2C0C1C2C2C0C6C7C0C0C0C0C0C4C3C1C0C7C5C0C1C5C1C5C1C8"
  mdtag <- getMD(a)
  expect_is(mdtag, "list")
  expect_length(mdtag$MDtag, 75)
  expect_length(mdtag$nucl.num, 75)
})

test_that("large insertion in MD", {
  a <- "2C0C3C7C2C7C1C0C1C1^ACGT33C1C2C2C"
  mdtag <- getMD(a)
  expect_is(mdtag, "list")
  expect_length(mdtag$MDtag, 31)
  expect_length(mdtag$nucl.num, 31)
})

test_that("single insertion in MD", {
  a <- "4C5C0C2C0C4C0G2C0C0C0C19C4^C2C2C0C1C4C5C10C5C5C4"
  mdtag <- getMD(a)
  expect_is(mdtag, "list")
  expect_length(mdtag$MDtag, 45)
  expect_length(mdtag$nucl.num, 45)
})

