context("test-bigwigs")

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

chrom.sizes <- get_data_path("hg19.chrom.sizes.mod")

sample_list <- try(readRDS(get_data_path("splitbams_output.rds")))
derASM <- calc_derivedasm(sample_list, cores = 1, verbose = T)

test_that("end to end make bigwigs", {
  #for one sample
  make_bigwig(colnames(derASM)[1], derASM, folder = ".", chromsizes.file = chrom.sizes)
  bw <- try(rtracklayer::import("der.ASM.CRC1.bw"))
  expect_s4_class(bw, "GRanges")
  
  #cleanup
  if(file.exists("der.ASM.CRC1.bedgraph")) file.remove("der.ASM.CRC1.sorted.bedgraph")
  if(file.exists("der.ASM.CRC1.bw")) file.remove("der.ASM.CRC1.bw")
  if(file.exists("der.ASM.CRC1.sorted.bedgraph")) file.remove("der.ASM.CRC1.sorted.bedgraph")
})
