context("split_bams")

DATA_PATH_DIR <- "../../inst/extdata/"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

#bam_files <- sapply(c("NORM1_chr19.bam", "CRC1_chr19.bam"), get_data_path)
bam_files <- get_data_path("NORM1_chr19_trim.bam")

#vcf_files <- sapply(c("NORM1.chr19.moretrim.vcf", "CRC1.chr19.moretrim.vcf"), get_data_path)
vcf_files <- get_data_path("NORM1.chr19.moretrim.vcf")

#sample_names <- c("NORM1", "CRC1")
sample_names <- "NORM1"

reference_file <- get_data_path("19.fa")
output_file <- "snp.table.NORM1.rds" #TODO make as an argument

split_bams(bam_files, vcf_files, sample_names, reference_file)

test_that("end to end split_bams", {
  split_bams(bam_files, vcf_files, sample_names, reference_file)
  GRanges_list <- try(readRDS(output_file))
  expect_type(GRanges_list, "list")

  #cleanup
  if(file.exists(output_file)) file.remove(output_file) #make function?
})

#expect_type(XX, "list")
#expect_s4_class(XX[[1]], "GRanges")
