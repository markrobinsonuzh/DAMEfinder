context("test-mdsplot")

DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

rds_files <- list.files(DATA_PATH_DIR, ".rds$")
rds_files <- get_data_path(rds_files)
derASM <- calc_derivedasm(rds_files, cores = 1, verbose = T)

test_that("end to end mdsplot", {
  p <- methyl_MDS_plot(derASM, color = c(rep("CRC",4), rep("NORM", 4)))
  expect_true(is.ggplot(p))
})
