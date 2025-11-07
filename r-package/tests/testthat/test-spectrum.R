test_that("spectrum post-hoc indicators can be output", {
  parameters <- read_parameters(test_path("testdata/spectrum_params.h5"))
  out <- run_model(parameters, "Spectrum")

  expected_output <- c(
    "p_deaths_nonaids_artpop", "p_deaths_nonaids_hivpop",
    "p_excess_deaths_nonaids_no_art", "p_excess_deaths_nonaids_on_art"
  )
  expect_true(all(expected_output %in% names(out)))
  # Always 0 for children
  expect_true(all(out$p_deaths_nonaids_artpop[1:15, , ] == 0))
  expect_true(all(out$p_deaths_nonaids_hivpop[1:15, , ] == 0))
  expect_true(any(out$p_deaths_nonaids_artpop > 0))
  expect_true(any(out$p_deaths_nonaids_hivpop > 0))

  expect_true(all(out$p_excess_deaths_nonaids_no_art[1:15, , ] == 0))
  expect_true(all(out$p_excess_deaths_nonaids_on_art[1:15, , ] == 0))
  expect_true(any(out$p_excess_deaths_nonaids_no_art > 0))
  expect_true(any(out$p_excess_deaths_nonaids_on_art > 0))
})

test_that("bwa test data returns 0 non-aids excess deaths", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  out <- run_model(parameters, "Spectrum")

  expected_output <- c(
    "p_deaths_nonaids_artpop", "p_deaths_nonaids_hivpop",
    "p_excess_deaths_nonaids_no_art", "p_excess_deaths_nonaids_on_art"
  )
  expect_true(all(expected_output %in% names(out)))

  # Always 0 with BWA test data
  expect_true(all(out$p_deaths_nonaids_artpop == 0))
  expect_true(all(out$p_deaths_nonaids_hivpop[1:15, , ] == 0))
  expect_true(any(out$p_deaths_nonaids_hivpop > 0))

  expect_true(all(out$p_excess_deaths_nonaids_no_art == 0))
  expect_true(all(out$p_excess_deaths_nonaids_on_art == 0))
})
