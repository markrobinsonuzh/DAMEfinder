context("test-methyl_circle_plot")
library(GenomicRanges)

DATA_PATH_DIR <- "../../inst/extdata/"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)


bam_files <- get_data_path("NORM1_chr19_trim.bam")

vcf_files <- get_data_path("NORM1.chr19.trim.vcf")

sample_names <- "NORM1"

reference_file <- get_data_path("19.fa")

snp <- GRanges(19, IRanges(388065, width = 1))


test_that("end to end methyl plot", {
  
  p <- methyl_circle_plot(snp = snp,
                          vcf.file = vcf_files,
                          bam.file = bam_files,
                          ref.file = reference_file,
                          letter.size = 3,
                          #dame = dame,
                          sample.name = sample_names,
                          point.size = 2)
  expect_true(is.ggplot(p))
  
  
  #cleanup
  if(file.exists("NORM1.RData")) file.remove("NORM1.RData")
})

test_that("add CpG site", {
  cpgsite <- GRanges(19, IRanges(387982, width = 1))
  p <- methyl_circle_plot(snp = snp,
                          vcf.file = vcf_files,
                          bam.file = bam_files,
                          ref.file = reference_file,
                          letter.size = 3,
                          #dame = dame,
                          sample.name = sample_names,
                          point.size = 2,
                          cpgsite = cpgsite)
  expect_true(is.ggplot(p))
  
  
  #cleanup
  if(file.exists("NORM1.RData")) file.remove("NORM1.RData")
})

test_that("end to end methyl_circle_plotCpG", {
  methyl_circle_plotCpG(cpgsite = cpgsite,
                        bam.file = bam_files,
                        ref.file = reference_file,
                        point.size = 2)
  expect_true(is.ggplot(p))
})

