context("split_bams")
#library(DAMEfinder)

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)


bam_files <- get_data_path("NORM1_chr19_trim.bam")

vcf_files <- get_data_path("NORM1.chr19.trim.vcf")
  
sample_names <- "NORM1"

reference_file <- get_data_path("19.fa")

test_that("end to end split_bams", {
  GRanges_list <- extract_bams(bam_files, vcf_files, sample_names, reference_file)
  expect_type(GRanges_list, "list")
})


test_that("output is GRangesList", {
  GRanges_list <- extract_bams(bam_files, vcf_files, sample_names, reference_file)
  expect_s4_class(GRanges_list$NORM1[[1]], "GRanges")
})



