context("test-methyl_circle_plot")
library(GenomicRanges)

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)


bam_files <- get_data_path("NORM1_chr19_trim.bam")

vcf_files <- get_data_path("NORM1.chr19.trim.vcf")

sample_names <- "NORM1"

#reference_file <- get_data_path("19.fa")
#Get reference file 
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
seqnames(genome) <- gsub("chr","",seqnames(genome))
reference_file <- DNAStringSet(genome[[19]], use.names = TRUE)
names(reference_file) <- 19

snp <- GRanges(19, IRanges(388065, width = 1))
cpgsite <- GRanges(19, IRanges(387982, width = 1))


test_that("end to end methyl plot", {
  
  p <- methyl_circle_plot(snp = snp,
                          vcfFile = vcf_files,
                          bamFile = bam_files,
                          refFile = reference_file,
                          letterSize = 3,
                          #dame = dame,
                          sampleName = sample_names,
                          pointSize = 2)
  expect_true(is.ggplot(p))
  
  
  #cleanup
  if(file.exists("NORM1.RData")) file.remove("NORM1.RData")
})

test_that("add CpG site", {
  p <- methyl_circle_plot(snp = snp,
                          vcfFile = vcf_files,
                          bamFile = bam_files,
                          refFile = reference_file,
                          letterSize = 3,
                          #dame = dame,
                          sampleName = sample_names,
                          pointSize = 2,
                          cpgsite = cpgsite)
  expect_true(is.ggplot(p))
  
  
  #cleanup
  if(file.exists("NORM1.RData")) file.remove("NORM1.RData")
})

test_that("end to end methyl_circle_plotCpG", {
  p <- methyl_circle_plotCpG(cpgsite = cpgsite,
                        bamFile = bam_files,
                        refFile = reference_file,
                        pointSize = 2)
  expect_true(is.ggplot(p))
})

