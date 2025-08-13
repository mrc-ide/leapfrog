test_that("initial state set up works as expected", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))

  out <- run_model(parameters, "HivFullAgeStratification", 1970L)

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_background_deaths", "p_hiv_pop",
      "p_hiv_pop_background_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 1))
  expect_equal(out$p_total_pop[, , 1], parameters$basepop[, , 1],
    ignore_attr = TRUE
  )

  ## Base year calculations will set values for
  ## p_total_pop, births and p_total_pop_background_deaths
  expect_true(out$births[1] != 0)
  expect_true(all(out$p_total_pop_background_deaths != 0))
  expect_equal(out$p_hiv_pop, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
  expect_equal(out$p_hiv_pop_background_deaths, array(rep(0, 81 * 2),
    dim = c(81, 2, 1)
  ))
  expect_equal(
    out$h_hiv_adult,
    array(rep(0, 7 * 66 * 2), dim = c(7, 66, 2, 1))
  )
  expect_equal(
    out$h_art_adult,
    array(rep(0, 3 * 7 * 66 * 2), dim = c(3, 7, 66, 2, 1))
  )
  expect_equal(
    out$h_hiv_deaths_no_art,
    array(rep(0, 7 * 66 * 2), dim = c(7, 66, 2, 1))
  )
  expect_equal(out$p_infections, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
  expect_equal(
    out$h_hiv_deaths_art,
    array(rep(0, 3 * 7 * 66 * 2), dim = c(3, 7, 66, 2, 1))
  )
  expect_equal(
    out$h_art_initiation,
    array(rep(0, 7 * 66 * 2), dim = c(7, 66, 2, 1))
  )
  expect_equal(out$p_hiv_deaths, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
})

test_that("initial state set up with coarse stratified HIV works as expected", {
  parameters <- read_parameters(test_path("testdata/adult_parms_coarse.h5"))

  out <- run_model(parameters, "HivCoarseAgeStratification", 1970L)

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_background_deaths", "p_hiv_pop",
      "p_hiv_pop_background_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 1))
  expect_equal(out$p_total_pop[, , 1], parameters$basepop[, , 1],
               ignore_attr = TRUE)

  ## Base year calculations will set values for
  ## p_total_pop, births and p_total_pop_background_deaths
  expect_true(out$births[1] != 0)
  expect_true(all(out$p_total_pop_background_deaths != 0))
  expect_equal(out$p_hiv_pop, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
  expect_equal(out$p_hiv_pop_background_deaths, array(rep(0, 81 * 2),
    dim = c(81, 2, 1)
  ))
  expect_equal(
    out$h_hiv_adult,
    array(rep(0, 7 * 9 * 2), dim = c(7, 9, 2, 1))
  )
  expect_equal(
    out$h_art_adult,
    array(rep(0, 3 * 7 * 9 * 2), dim = c(3, 7, 9, 2, 1))
  )
  expect_equal(
    out$h_hiv_deaths_no_art,
    array(rep(0, 7 * 9 * 2), dim = c(7, 9, 2, 1))
  )
  expect_equal(out$p_infections, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
  expect_equal(
    out$h_hiv_deaths_art,
    array(rep(0, 3 * 7 * 9 * 2), dim = c(3, 7, 9, 2, 1))
  )
  expect_equal(
    out$h_art_initiation,
    array(rep(0, 7 * 9 * 2), dim = c(7, 9, 2, 1))
  )
  expect_equal(out$p_hiv_deaths, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
})

test_that("model for 1 time step has looped", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))

  out <- run_model(parameters, "HivFullAgeStratification", 1971)

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_background_deaths", "p_hiv_pop",
      "p_hiv_pop_background_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 1))
  expect_true(all(out$p_total_pop > 100))
  expect_true(out$births > 0) ## a simulation has been run this is not still 0
  expect_equal(dim(out$p_total_pop_background_deaths), c(81, 2, 1))
  expect_true(all(out$p_total_pop_background_deaths > 0))
  ## These are going to stay 0 as no p_infections after just 1 year has been run
  expect_true(all(out$p_hiv_pop == 0))
  expect_true(all(out$p_hiv_pop_background_deaths == 0))
  expect_true(all(out$h_hiv_adult == 0))
  expect_true(all(out$h_art_adult == 0))
  expect_true(all(out$h_hiv_deaths_no_art == 0))
  expect_true(all(out$p_infections == 0))
  expect_true(all(out$h_hiv_deaths_art == 0))
  expect_true(all(out$h_art_initiation == 0))
  expect_true(all(out$p_hiv_deaths == 0))
})

test_that("model can be run for all years", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))

  out <- run_model(parameters)

  ## No HIV population < age 15
  expect_true(all(out$p_hiv_pop[1:15, , ] < 1e-20))
  expect_true(all(out$p_hiv_pop[1:15, , ] > -1e-20))
  expect_true(all(out$p_hiv_pop_background_deaths[1:16, , ] == 0))
  expect_true(all(out$p_infections[1:15, , ] == 0))

  ## There is HIV population after age 15
  expect_true(all(out$p_hiv_pop[16:nrow(out$p_hiv_pop), , 61] > 0))
  ## Natural deaths start at index 17 as no deaths in first HIV population
  ## projection as they are calculated from
  ## the no of HIV +ve in previous year - is this right?
  expect_true(
    all(out$p_hiv_pop_background_deaths[17:nrow(out$p_hiv_pop), , 61] != 0)
  )
  ## Some of older ages can be 0 p_infections, so check the middle chunk
  expect_true(all(out$p_infections[16:70, , 61] > 0))

  expect_true(all(out$h_hiv_adult[, , , 61] != 0))
  expect_true(all(out$h_art_adult[, , , , 61] != 0))
  expect_true(all(out$h_art_initiation[, , , 61] != 0))

  ## Outputs cannot be negative
  expect_true(all(out$p_total_pop >= 0))
  expect_true(all(out$births >= 0))
  expect_true(all(out$p_total_pop_background_deaths >= 0))
  expect_true(all(out$p_hiv_pop >= 0))
  expect_true(all(out$p_hiv_pop_background_deaths >= 0))
  expect_true(all(out$h_hiv_adult >= 0))
  expect_true(all(out$h_art_adult >= 0))
  expect_true(all(out$h_hiv_deaths_no_art >= 0))
  expect_true(all(out$p_infections >= 0))
  expect_true(all(out$h_art_initiation >= 0))
  expect_true(all(out$h_hiv_deaths_art >= 0))
  expect_true(all(out$p_hiv_deaths >= 0))
})

test_that("model can be run with ART initiation", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))
  ## Set time ART start to some value lower than no of years in projection
  parameters[["t_ART_start"]] <- 20L

  expect_silent(run_model(parameters))
})

test_that("model can be run twice on the same data", {
  ## Regression test as we saw the 2nd run failing as the first fit
  ## was modifying the R stored data causing the 2nd run on the same
  ## data to read from an index of -1
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))

  out <- run_model(parameters)
  out2 <- run_model(parameters)
  expect_identical(out, out2)
})

test_that("child model can be run twice on the same data", {
  ## Regression test as we saw the 2nd run failing as the first fit
  ## was modifying the R stored data causing the 2nd run on the same
  ## data to read from an index of -1
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))

  out <- run_model(parameters, "ChildModel", 1970:2030)
  out2 <- run_model(parameters, "ChildModel", 1970:2030)
  expect_identical(out, out2)
})

test_that("error thrown if model run with invalid configuration", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))

  expect_error(
    run_model(parameters, "HivFineAgeStratification", 2030L),
    paste("Invalid configuration: 'HivFineAgeStratification'.",
          "It must be one of: 'DemographicProjection'")
  )
})

test_that("error thrown if size of stratified data does not match expected", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))
  parameters[["cd4_mort"]] <- rep(1, 3)

  expect_error(
    run_model(parameters, "HivFullAgeStratification", 2030L),
    "Invalid size of data for 'cd4_mort', expected 924 got 3"
  )
})

test_that("error if trying to save output from before projection start", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))

  expect_error(
    run_model(parameters, "ChildModel", 1965:2030),
    paste("Trying to output for year: '1965' which",
          "is before the projection start year: '1970'.")
  )
})

test_that("error thrown if invalid projection period set", {
  parameters <- read_parameters(test_path("testdata/adult_parms.h5"))
  parameters$projection_period <- "unknown"

  expect_error(
    run_model(parameters, "HivFullAgeStratification", 2030L),
    "Invalid projection period: 'unknown'. Allowed values are: 'midyear' or 'calendar'."
  )
})
