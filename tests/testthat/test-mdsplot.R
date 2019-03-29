context("test-mdsplot")

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

sample_list <- try(readRDS(get_data_path("splitbams_output.rds")))

derASM <- calc_derivedasm(sample_list, cores = 1, verbose = T)

test_that("end to end mdsplot", {
  p <- methyl_MDS_plot(derASM, color = c(rep("CRC",4), rep("NORM", 4)))
  expect_true(is.ggplot(p))
})
