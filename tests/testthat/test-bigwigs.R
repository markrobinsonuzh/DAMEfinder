context("test-bigwigs")

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

chrom.sizes <- get_data_path("hg19.chrom.sizes.mod")

rds_files <- list.files(DATA_PATH_DIR, ".rds$")
rds_files <- get_data_path(rds_files)
derASM <- calc_derivedasm(rds_files, cores = 1, verbose = T)

test_that("end to end make bigwigs", {
  make_bigwig(colnames(derASM)[1], derASM, folder = ".", chromsizes.file = chrom.sizes)
  bw <- try(rtracklayer::import(get_data_path("der.ASM.CRC1.bw")))
  expect_type(bw, "GRanges")
  
  #cleanup
  if(file.exists("der.ASM.CRC1.bedgraph")) file.remove("der.ASM.CRC1.sorted.bedgraph")
  if(file.exists("der.ASM.CRC1.bw")) file.remove("der.ASM.CRC1.bw")
  if(file.exists("der.ASM.CRC1.sorted.bedgraph")) file.remove("der.ASM.CRC1.sorted.bedgraph")
})
