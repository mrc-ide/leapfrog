test_that("demographic projection can be run", {
  parameters <- read_parameters(test_path("testdata/adult_parms_full.h5"))

  out <- run_model(parameters, "DemographicProjection", 1970:2030)

  expect_setequal(
    names(out),
    c(
      "p_totpop", "births", "p_background_deaths_totpop"
    )
  )
  expect_equal(dim(out$p_totpop), c(81, 2, 61))
  expect_true(all(out$p_totpop > 100))
  expect_true(out$births[61] > 0) ## a simulation has been run this is not still 0
  expect_equal(dim(out$p_background_deaths_totpop), c(81, 2, 61))
  expect_true(all(out$p_background_deaths_totpop[, , 61] > 0))
})
