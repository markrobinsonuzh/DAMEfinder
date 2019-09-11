context("test-mdsplot")

data(extractbams_output)

derASM <- calc_derivedasm(extractbams_output, cores = 1, verbose = TRUE)

test_that("end to end mdsplot", {
  p <- methyl_MDS_plot(derASM, group = c(rep("CRC",4), rep("NORM", 4)))
  expect_true(is.ggplot(p))
})
