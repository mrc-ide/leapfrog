test_that("initial state set up works as expected", {
  parameters <- read_parameters(test_path("testdata/adult_parms_coarse.h5"))

  parameters$rvec <- rep(0.1, parameters$hts_per_year * 61)

  out <- run_model(parameters, "HivCoarseAgeFit")

  # Just a smoke test for now, need to flesh it out
  expect_true(any(out$p_totpop != 0.0))
})
