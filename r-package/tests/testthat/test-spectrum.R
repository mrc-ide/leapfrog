test_that("spectrum post-hoc indicators can be output", {
  parameters <- read_parameters(test_path("testdata/spectrum_params.h5"))

  out <- run_model(parameters, "Spectrum")

  expect_true(all(c("p_deaths_nonaids_artpop", "p_deaths_nonaids_hivpop") %in% names(out)))
  # Always 0 for children
  expect_true(all(out$p_deaths_nonaids_artpop[1:15, , ] == 0))
  expect_true(all(out$p_deaths_nonaids_hivpop[1:15, , ] == 0))
  expect_true(any(out$p_deaths_nonaids_artpop > 0))
  expect_true(any(out$p_deaths_nonaids_hivpop > 0))
})
