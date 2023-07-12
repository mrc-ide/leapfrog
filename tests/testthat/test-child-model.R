test_that("child model can be run for all years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_child.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_child.rds"))

  expect_silent(out <- run_model(demp, parameters, NULL, NULL, 0:60))

  expect_true(all(out$hc_hiv_pop[, , , , 1:which(1970:2030 == 1979)] == 0))

  ## All 10 as seeded with 100 infections which are distributed over 5 age
  ## groups and genders evenly
  expect_true(all(out$infections[1:5, , which(1970:2030 == 1980)] == 10))
  expect_true(all(out$hiv_population[1:5, , which(1970:2030 == 1980)] == 10))
  expect_true(all(out$hc_hiv_pop[1, 1, 1:5, , which(1970:2030 == 1980)] == 10))

  ## Only seeding nosocomial infections in year 1980, infections every
  ## other year should be 0
  expect_true(all(out$infections[1:15, , which(1970:2030 != 1980)] == 0))
})
