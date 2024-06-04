test_that("initial state set up works as expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  out <- run_model(demp, parameters, 1970L, 0L, run_child_model = FALSE)

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths", "p_hiv_pop",
      "p_hiv_pop_natural_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 1))
  expect_equal(out$p_total_pop[, , 1], demp$basepop[, , 1],
    ignore_attr = TRUE
  )

  expect_equal(out$births[1], 0)
  expect_equal(out$p_total_pop_natural_deaths, array(rep(0, 81 * 2),
    dim = c(81, 2, 1)
  ))
  expect_equal(out$p_hiv_pop, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
  expect_equal(out$p_hiv_pop_natural_deaths, array(rep(0, 81 * 2),
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
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  out <- run_model(demp, parameters, 1970L, 0L,
    hiv_age_stratification = "coarse",
    run_child_model = FALSE
  )

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths", "p_hiv_pop",
      "p_hiv_pop_natural_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 1))
  expect_equal(out$p_total_pop[, , 1], demp$basepop[, , 1], ignore_attr = TRUE)

  expect_equal(out$births[1], 0)
  expect_equal(out$p_total_pop_natural_deaths, array(rep(0, 81 * 2),
    dim = c(81, 2, 1)
  ))
  expect_equal(out$p_hiv_pop, array(rep(0, 81 * 2), dim = c(81, 2, 1)))
  expect_equal(out$p_hiv_pop_natural_deaths, array(rep(0, 81 * 2),
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
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  out <- run_model(demp, parameters, 1970:1971, 10L, 1971L,
    run_child_model = FALSE
  )

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths", "p_hiv_pop",
      "p_hiv_pop_natural_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths"
    )
  )
  expect_equal(dim(out$p_total_pop), c(81, 2, 1))
  expect_true(all(out$p_total_pop > 100))
  expect_true(out$births > 0) ## a simulation has been run this is not still 0
  expect_equal(dim(out$p_total_pop_natural_deaths), c(81, 2, 1))
  expect_true(all(out$p_total_pop_natural_deaths > 0))
  ## These are going to stay 0 as no p_infections after just 1 year has been run
  expect_true(all(out$p_hiv_pop == 0))
  expect_true(all(out$p_hiv_pop_natural_deaths == 0))
  expect_true(all(out$h_hiv_adult == 0))
  expect_true(all(out$h_art_adult == 0))
  expect_true(all(out$h_hiv_deaths_no_art == 0))
  expect_true(all(out$p_infections == 0))
  expect_true(all(out$h_hiv_deaths_art == 0))
  expect_true(all(out$h_art_initiation == 0))
  expect_true(all(out$p_hiv_deaths == 0))
})

test_that("model can be run for all years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  out <- run_model(demp, parameters, NULL, NULL,
    run_child_model = FALSE
  )

  ## No HIV population < age 15
  expect_true(all(out$p_hiv_pop[1:15, , ] < 1e-20))
  expect_true(all(out$p_hiv_pop[1:15, , ] > -1e-20))
  expect_true(all(out$p_hiv_pop_natural_deaths[1:16, , ] == 0))
  expect_true(all(out$p_infections[1:15, , ] == 0))

  ## There is HIV population after age 15
  expect_true(all(out$p_hiv_pop[16:nrow(out$p_hiv_pop), , 61] > 0))
  ## Natural deaths start at index 17 as no deaths in first HIV population
  ## projection as they are calculated from
  ## the no of HIV +ve in previous year - is this right?
  expect_true(
    all(out$p_hiv_pop_natural_deaths[17:nrow(out$p_hiv_pop), , 61] != 0)
  )
  ## Some of older ages can be 0 p_infections, so check the middle chunk
  expect_true(all(out$p_infections[16:70, , 61] > 0))

  expect_true(all(out$h_hiv_adult[, , , 61] != 0))
  expect_true(all(out$h_art_adult[, , , , 61] != 0))
  expect_true(all(out$h_art_initiation[, , , 61] != 0))

  ## Outputs cannot be negative
  expect_true(all(out$p_total_pop >= 0))
  expect_true(all(out$births >= 0))
  expect_true(all(out$p_total_pop_natural_deaths >= 0))
  expect_true(all(out$p_hiv_pop >= 0))
  expect_true(all(out$p_hiv_pop_natural_deaths >= 0))
  expect_true(all(out$h_hiv_adult >= 0))
  expect_true(all(out$h_art_adult >= 0))
  expect_true(all(out$h_hiv_deaths_no_art >= 0))
  expect_true(all(out$p_infections >= 0))
  expect_true(all(out$h_art_initiation >= 0))
  expect_true(all(out$p_infections >= 0))
  expect_true(all(out$p_infections >= 0))
})

test_that("model can be run with ART initiation", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))
  ## Set time ART start to some value lower than no of years in projection
  parameters[["t_ART_start"]] <- 20L

  expect_silent(run_model(demp, parameters, NULL, NULL,
    run_child_model = FALSE
  ))
})

test_that("model can be run twice on the same data", {
  ## Regression test as we saw the 2nd run failing as the first fit
  ## was modifying the R stored data causing the 2nd run on the same
  ## data to read from an index of -1
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  out <- run_model(demp, parameters, NULL, NULL,
    run_child_model = FALSE
  )
  out2 <- run_model(demp, parameters, NULL, NULL,
    run_child_model = FALSE
  )
  expect_identical(out, out2)
})

test_that("child model can be run twice on the same data", {
  ## Regression test as we saw the 2nd run failing as the first fit
  ## was modifying the R stored data causing the 2nd run on the same
  ## data to read from an index of -1
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")

  out <- run_model(input$demp, input$parameters, NULL, NULL,
    run_child_model = TRUE
  )
  out2 <- run_model(input$demp, input$parameters, NULL, NULL,
    run_child_model = TRUE
  )
  expect_identical(out, out2)
})

test_that("error thrown if trying to run model for more than max years", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  expect_error(
    run_model(demp, parameters, 1970:2050, 10L, 1970:2050),
    "No of years > max years of 61"
  )
})

test_that("error thrown if model run with invalid HIV stratification", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  expect_error(
    run_model(demp, parameters, NULL, NULL, 2030L,
      hiv_age_stratification = "fine"
    ),
    "hiv_age_stratification must be one of 'full', 'coarse', got 'fine'"
  )
})

test_that("error thrown if size of stratified data does not match expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))
  parameters[["cd4_mort_full"]] <- rep(1, 3)

  expect_error(
    run_model(demp, parameters, NULL, NULL, 2030L,
      hiv_age_stratification = "full"
    ),
    "Invalid size of data for 'cd4_mort', expected 924 got 3"
  )
})

test_that("error thrown if trying to save output from invalid steps", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")

  expect_error(
    run_model(input$demp, input$parameters, NULL, NULL, -1L),
    "Invalid output step '-1'. Can only output one of the simulation years."
  )

  expect_error(
    run_model(input$demp, input$parameters, 1970:1980, NULL, 1981:1982),
    paste(
      "Invalid output steps '1981', '1982'.",
      "Can only output one of the simulation years."
    )
  )
})
