test_that("initial state set up works as expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  out <- run_base_model(demp, parameters, 0L, 0L)

  expect_setequal(
    names(out),
    c(
      "total_population", "births", "natural_deaths", "hiv_population",
      "hiv_natural_deaths", "hiv_strat_adult", "art_strat_adult",
      "aids_deaths_no_art", "infections", "aids_deaths_art"
    )
  )
  expect_equal(dim(out$total_population), c(81, 2))
  expect_equal(out$total_population, demp$basepop[, , 1], ignore_attr = TRUE)

  expect_equal(out$births, 0)
  expect_equal(out$natural_deaths, matrix(0, nrow = 81, ncol = 2))
  expect_equal(out$hiv_population, matrix(0, nrow = 81, ncol = 2))
  expect_equal(out$hiv_natural_deaths, matrix(0, nrow = 81, ncol = 2))
  expect_equal(
    out$hiv_strat_adult,
    array(rep(0, 7 * 66 * 2), dim = c(7, 66, 2))
  )
  expect_equal(
    out$art_strat_adult,
    array(rep(0, 3 * 7 * 66 * 2), dim = c(3, 7, 66, 2))
  )
  expect_equal(
    out$aids_deaths_no_art,
    array(rep(0, 7 * 66 * 2), dim = c(7, 66, 2))
  )
  expect_equal(out$infections, matrix(0, nrow = 81, ncol = 2))
  expect_equal(
    out$aids_deaths_art,
    array(rep(0, 3 * 7 * 66 * 2), dim = c(3, 7, 66, 2))
  )
})

test_that("initial state set up with coarse stratified HIV works as expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  out <- run_base_model(demp, parameters, 0L, 0L, hiv_age_stratification = "coarse")

  expect_setequal(
    names(out),
    c(
      "total_population", "births", "natural_deaths", "hiv_population",
      "hiv_natural_deaths", "hiv_strat_adult", "art_strat_adult",
      "aids_deaths_no_art", "infections", "aids_deaths_art"
    )
  )
  expect_equal(dim(out$total_population), c(81, 2))
  expect_equal(out$total_population, demp$basepop[, , 1], ignore_attr = TRUE)

  expect_equal(out$births, 0)
  expect_equal(out$natural_deaths, matrix(0, nrow = 81, ncol = 2))
  expect_equal(out$hiv_population, matrix(0, nrow = 81, ncol = 2))
  expect_equal(out$hiv_natural_deaths, matrix(0, nrow = 81, ncol = 2))
  expect_equal(
    out$hiv_strat_adult,
    array(rep(0, 7 * 9 * 2), dim = c(7, 9, 2))
  )
  expect_equal(
    out$art_strat_adult,
    array(rep(0, 3 * 7 * 9 * 2), dim = c(3, 7, 9, 2))
  )
  expect_equal(
    out$aids_deaths_no_art,
    array(rep(0, 7 * 9 * 2), dim = c(7, 9, 2))
  )
  expect_equal(out$infections, matrix(0, nrow = 81, ncol = 2))
  expect_equal(
    out$aids_deaths_art,
    array(rep(0, 3 * 7 * 9 * 2), dim = c(3, 7, 9, 2))
  )
})

test_that("model for 1 time step has looped", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  out <- run_base_model(demp, parameters, 1L, 10L)

  expect_setequal(
    names(out),
    c(
      "total_population", "births", "natural_deaths", "hiv_population",
      "hiv_natural_deaths", "hiv_strat_adult", "art_strat_adult",
      "aids_deaths_no_art", "infections", "aids_deaths_art"
    )
  )
  expect_equal(dim(out$total_population), c(81, 2))
  expect_true(all(out$total_population > 100))
  expect_true(out$births > 0) ## a simulation has been run this is not still 0
  expect_equal(dim(out$natural_deaths), c(81, 2))
  expect_true(all(out$natural_deaths > 0))
  ## These are going to stay 0 as no infections after just 1 year has been run
  expect_true(all(out$hiv_population == 0))
  expect_true(all(out$hiv_natural_deaths == 0))
  expect_true(all(out$hiv_strat_adult == 0))
  expect_true(all(out$art_strat_adult == 0))
  expect_true(all(out$aids_deaths_no_art == 0))
  expect_true(all(out$infections == 0))
  ## Machine precision can mean these might not come out as nice round 0s
  expect_true(all(out$aids_deaths_art < 1e-10))
  expect_true(all(out$aids_deaths_art > -1e-10))
})

test_that("model can be run for all years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  expect_error(out <- run_base_model(demp, parameters, NULL, NULL), NA)

  ## No HIV population < age 15
  expect_true(all(out$hiv_population[1:15, ] == 0))
  expect_true(all(out$hiv_natural_deaths[1:16, ] == 0))
  expect_true(all(out$infections[1:15, ] == 0))

  ## There is HIV population after age 15
  expect_true(all(out$hiv_population[16:nrow(out$hiv_population), ] > 0))
  ## Natural deaths start at index 17 as no deaths in first HIV population projection as they are calculated from
  ## the no of HIV +ve in previous year - is this right?
  expect_true(all(out$hiv_natural_deaths[17:nrow(out$hiv_population), ] != 0))
  ## Some of older ages can be 0 infections, so check the middle chunk
  expect_true(all(out$infections[16:70, ] != 0))

  ## HIV and ART strat still 0, we are not adding to these yet
  expect_true(all(out$hiv_strat_adult == 0))
  expect_true(all(out$art_strat_adult == 0))
  expect_equal(
    out$aids_deaths_no_art,
    array(rep(0, 7 * 66 * 2), dim = c(7, 66, 2))
  )
})

test_that("model can be run with ART initiation", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))
  ## Set time ART start to some value lower than no of years in projection
  parameters[["t_ART_start"]] <- 20L

  expect_error(run_base_model(demp, parameters, NULL, NULL), NA)
})

test_that("error thrown if trying to run model for more than max years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  expect_error(
    run_base_model(demp, parameters, 100L, 10L),
    "No of years > max years of 60"
  )
})

test_that("error thrown if model run with invalid HIV stratification", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  expect_error(
    run_base_model(demp, parameters, NULL, NULL, hiv_age_stratification = "fine"),
    "Invalid HIV age stratification must be 'full' or 'coarse' got 'fine'."
  )
})

test_that("error thrown if size of stratified data does not match expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))
  parameters[["hAG_SPAN_full"]] <- rep(1, 3)

  expect_error(
    run_base_model(demp, parameters, NULL, NULL, hiv_age_stratification = "full"),
    "Invalid size of data, expected 66 got 3"
  )
})
