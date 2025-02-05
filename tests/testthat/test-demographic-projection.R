test_that("demographic projection can be run", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  out <- run_model(demp, parameters, NULL, NULL,
    run_hiv_simulation = FALSE,
    run_child_model = FALSE
  )

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
