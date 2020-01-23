context("test-mdsplot")

data(readtuples_output)
ASM <- calc_asm(readtuples_output)

test_that("end to end mdsplot", {
  p <- methyl_MDS_plot(ASM, group = c(rep("CRC",3), rep("NORM", 2)))
  expect_true(is.ggplot(p))
})
