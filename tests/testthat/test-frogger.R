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

test_that("model can be run for all years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, NULL), NA)
})

test_that("error thrown if trying to run model for more than max years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, 100L), NA)
})

test_that("try only 1 gender in input data?", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, 100L), NA)
})

test_that("try no migration (does it work)", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, 100L), NA)
})

test_that("different births or survival rate?", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))

  expect_error(run_base_model(demp, 100L), NA)
})