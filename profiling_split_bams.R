#install.packages("profvis")
library(profvis)

DATA_PATH_DIR <- "inst/extdata"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

#bam_files <- sapply(c("NORM1_chr19.bam", "CRC1_chr19.bam"), get_data_path)
bam_files <- get_data_path("NORM1_chr19_trim.bam")

#vcf_files <- sapply(c("NORM1.chr19.moretrim.vcf", "CRC1.chr19.moretrim.vcf"), get_data_path)
vcf_files <- get_data_path("NORM1.chr19.moretrim.vcf")

#sample_names <- c("NORM1", "CRC1")
sample_names <- "NORM1"

reference_file <- get_data_path("19.fa")
output_file <- "snp.table.NORM1.rds"

profvis({
  split_bams(bam_files, vcf_files, sample_names, reference_file)
})

# split_bams(bam_files, vcf_files, sample_names, reference_file)
# open(FaFile(reference_file))
# 
# fi <- FaFile(reference_file)
# 
# fl <- system.file("extdata", "ce2dict1.fa", package="Rsamtools", mustWork=TRUE)
# fa <- open(FaFile(fl))  
