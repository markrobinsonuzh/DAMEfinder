context("split_bams")
#This is a long test

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)


bam_files <- get_data_path("NORM1_chr19_trim.bam")

vcf_files <- get_data_path("NORM1.chr19.trim.vcf")

sample_names <- "NORM1"

#This fasta is from ftp://ftp.ensembl.org/pub/grch37/release-91/fasta/homo_sapiens/dna/
#reference_file <- get_data_path("19.fa")

#Get reference file 
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
seqnames(genome) <- gsub("chr","",seqnames(genome))
dna <- DNAStringSet(genome[[19]], use.names = TRUE)
names(dna) <- 19


test_that("output is GRangesList", {
  GRanges_list <- extract_bams(bam_files, vcf_files, sample_names, dna)
  expect_s4_class(GRanges_list$NORM1[[1]], "GRanges")
})



