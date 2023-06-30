test_that("demographic model is correct", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))

  actual <- run_model(demp, parameters, NULL, 0L, 0:60)

  expected <- readRDS(test_path("testdata/fit_demography.rds"))

  expect_equal(actual$total_population, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(as.numeric(actual$births), expected$births)
  expect_equal(actual$natural_deaths, expected$natdeaths)
  expect_equal(actual$hiv_population, expected$hivpop1)
  expect_equal(actual$hiv_natural_deaths, expected$natdeaths_hivpop)
  expect_equal(actual$hiv_strat_adult, expected$hivstrat_adult)
  expect_equal(actual$art_strat_adult, expected$artstrat_adult)
  expect_equal(actual$aids_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$infections, expected$infections)
  expect_equal(actual$aids_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$art_initiation, expected$artinit)
  expect_equal(actual$hiv_deaths, expected$hivdeaths)
})


test_that("model agrees with leapfrog impl", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))
  actual <- run_model(demp, parameters, NULL, NULL, 0:60)

  expected <- readRDS(test_path("testdata/leapfrog_fit.rds"))

  expect_equal(actual$total_population, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(as.numeric(actual$births), expected$births)
  expect_equal(actual$natural_deaths, expected$natdeaths)
  expect_equal(actual$hiv_population, expected$hivpop1)
  expect_equal(actual$hiv_natural_deaths, expected$natdeaths_hivpop)
  expect_equal(actual$hiv_strat_adult, expected$hivstrat_adult)
  expect_equal(actual$art_strat_adult, expected$artstrat_adult)
  expect_equal(actual$aids_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$infections, expected$infections)
  expect_equal(actual$aids_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$art_initiation, expected$artinit)
  expect_equal(actual$hiv_deaths, expected$hivdeaths)
})

test_that("model agrees with leapfrog impl", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))
  actual <- run_model(demp, parameters, NULL, NULL, 0:60, hiv_age_stratification = "coarse")

  expected <- readRDS(test_path("testdata/leapfrog_fit_coarse.rds"))

  expect_equal(actual$total_population, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(as.numeric(actual$births), expected$births)
  expect_equal(actual$natural_deaths, expected$natdeaths)
  expect_equal(actual$hiv_population, expected$hivpop1)
  expect_equal(actual$hiv_natural_deaths, expected$natdeaths_hivpop)
  expect_equal(actual$hiv_strat_adult, expected$hivstrat_adult)
  expect_equal(actual$art_strat_adult, expected$artstrat_adult)
  expect_equal(actual$aids_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$infections, expected$infections)
  expect_equal(actual$aids_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$art_initiation, expected$artinit)
  expect_equal(actual$hiv_deaths, expected$hivdeaths)
})
