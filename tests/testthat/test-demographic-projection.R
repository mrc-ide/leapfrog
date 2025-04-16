test_that("demographic projection can be run", {
  parameters <- readRDS(test_path("testdata/adult_parms.rds"))

  out <- run_model(parameters, "DemographicProjection", 1970:2030)

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 61))
  expect_true(all(out$p_total_pop > 100))
  expect_true(out$births[61] > 0) ## a simulation has been run this is not still 0
  expect_equal(dim(out$p_total_pop_natural_deaths), c(81, 2, 61))
  expect_true(all(out$p_total_pop_natural_deaths[, , 61] > 0))
})
