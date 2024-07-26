test_that("Invalid options throw error", {

  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_demproj-only_spectrum-v6.28_2023-12-12.PJNZ")
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)

  demp2$projection_period <- "bogus"
  expect_error(leapfrogR(demp2, hivp2),
               'projection_period "bogus" not found. Please select "midyear" or "calendar"')
})
