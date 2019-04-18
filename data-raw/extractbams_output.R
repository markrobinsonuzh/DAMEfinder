library(DAMEfinder)

DATA_PATH_DIR <- "inst/extdata"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

bam_files <- sapply(c("CRC1_chr19_trim.bam",
                      "moredata/CRC2_chr19_trim.bam",
                      "moredata/CRC3_chr19_trim.bam",
                      "moredata/CRC4_chr19_trim.bam",
                      "NORM1_chr19_trim.bam",
                      "moredata/NORM2_chr19_trim.bam",
                      "moredata/NORM3_chr19_trim.bam",
                      "moredata/NORM4_chr19_trim.bam"),get_data_path,
                    USE.NAMES = FALSE)

vcf_files <- sapply(c("CRC1.chr19.trim.vcf",
                      "moredata/CRC2.chr19.trim.vcf",
                      "moredata/CRC3.chr19.trim.vcf",
                      "moredata/CRC4.chr19.trim.vcf",
                      "NORM1.chr19.trim.vcf",
                      "moredata/NORM2.chr19.trim.vcf",
                      "moredata/NORM3.chr19.trim.vcf",
                      "moredata/NORM4.chr19.trim.vcf"),get_data_path,
                    USE.NAMES = FALSE)

sample_names <- c("CRC1","CRC2","CRC3","CRC4","NORM1","NORM2","NORM3","NORM4")

reference_file <- get_data_path("19.fa")

extractbams_output <- extract_bams(bam_files, vcf_files, sample_names, 
                                 reference_file)

usethis::use_data(extractbams_output, overwrite = TRUE, compress = 'xz')
