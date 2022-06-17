test_that("initial state set up works as expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  out <- run_base_model(demp, 0L)

  expect_setequal(names(out), c("total_population", "births", "natural_deaths"))
  expect_equal(dim(out$total_population), c(81, 2))
  expect_equal(out$total_population, demp$basepop[, , 1], ignore_attr = TRUE)

  expect_equal(out$births, 0)

  expect_equal(dim(out$natural_deaths), c(81, 2))
  expect_equal(out$natural_deaths, matrix(0, nrow = 81, ncol = 2))
})

test_that("model for 1 time step has not looped", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  out <- run_base_model(demp, 1L)

  expect_setequal(names(out), c("total_population", "births", "natural_deaths"))
  expect_equal(dim(out$total_population), c(81, 2))
  expect_true(all(out$total_population > 100))

  expect_true(out$births > 0) ## a simulation has been run this is not still 0

  expect_equal(dim(out$natural_deaths), c(81, 2))
  expect_true(all(out$natural_deaths > 0))
})

test_that("model can be run for all years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, NULL), NA)
})

test_that("error thrown if trying to run model for more than max years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, 100L), "No of years > max years of 60")
})
