test_that("demographic model is correct", {
  parameters <- readRDS(test_path("testdata/adult_parms.rds"))

  actual <- run_model(parameters, "DemographicProjection")

  expected <- readRDS(test_path("testdata/fit_demography.rds"))

  expect_equal(actual$p_total_pop, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(as.numeric(actual$births), expected$births)
  expect_equal(actual$p_total_pop_background_deaths, expected$natdeaths)
})


test_that("model agrees with leapfrog impl", {
  parameters <- readRDS(test_path("testdata/adult_parms.rds"))
  actual <- run_model(parameters)

  expected <- readRDS(test_path("testdata/leapfrog_fit.rds"))

  expect_equal(actual$p_total_pop, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(as.numeric(actual$births), expected$births)
  expect_equal(actual$p_total_pop_background_deaths, expected$natdeaths)
  expect_equal(actual$p_hiv_pop, expected$hivpop1)
  expect_equal(actual$p_hiv_pop_background_deaths, expected$natdeaths_hivpop)
  expect_equal(actual$h_hiv_adult, expected$hivstrat_adult)
  expect_equal(actual$h_art_adult, expected$artstrat_adult)
  expect_equal(actual$h_hiv_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$p_infections, expected$infections)
  expect_equal(actual$h_hiv_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$h_art_initiation, expected$artinit)
  expect_equal(actual$p_hiv_deaths, expected$hivdeaths)
})

test_that("model agrees with leapfrog impl", {
  parameters <- readRDS(test_path("testdata/adult_parms.rds"))
  actual <- run_model(parameters, "HivCoarseAgeStratification")

  expected <- readRDS(test_path("testdata/leapfrog_fit_coarse.rds"))

  expect_equal(actual$p_total_pop, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(as.numeric(actual$births), expected$births)
  expect_equal(actual$p_total_pop_background_deaths, expected$natdeaths)
  expect_equal(actual$p_hiv_pop, expected$hivpop1)
  expect_equal(actual$p_hiv_pop_background_deaths, expected$natdeaths_hivpop)
  expect_equal(actual$h_hiv_adult, expected$hivstrat_adult)
  expect_equal(actual$h_art_adult, expected$artstrat_adult)
  expect_equal(actual$h_hiv_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$p_infections, expected$infections)
  expect_equal(actual$h_hiv_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$h_art_initiation, expected$artinit)
  expect_equal(actual$p_hiv_deaths, expected$hivdeaths)
})
