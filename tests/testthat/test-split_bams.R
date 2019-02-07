context("split_bams")
#library(DAMEfinder)

DATA_PATH_DIR <- "../../inst/extdata/"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)


bam_files <- get_data_path("NORM1_chr19_trim.bam")

vcf_files <- get_data_path("NORM1.chr19.trim.vcf")
  
sample_names <- "NORM1"

reference_file <- get_data_path("19.fa")

output_file <- "snp.table.NORM1.rds"

test_that("end to end split_bams", {
  split_bams(bam_files, vcf_files, sample_names, reference_file)
  GRanges_list <- try(readRDS(output_file))
  expect_type(GRanges_list, "list")

  #cleanup
  if(file.exists(output_file)) file.remove(output_file) #make function?
})


test_that("output is GRangesList", {
  split_bams(bam_files, vcf_files, sample_names, reference_file)
  GRanges_list <- try(readRDS(output_file))
  expect_s4_class(GRanges_list[[1]], "GRanges")

  #cleanup
  if(file.exists(output_file)) file.remove(output_file)
})



