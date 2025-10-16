test_that("demographic model is correct", {
  parameters <- read_parameters(test_path("testdata/adult_parms_full.h5"))

  actual <- run_model(parameters, "DemographicProjection")

  expected <- read_parameters(test_path("testdata/fit_demography.h5"))

  expect_equal(actual$p_totpop, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(actual$births, expected$births)
  expect_equal(actual$p_background_deaths_totpop, expected$natdeaths)
})


test_that("model agrees with leapfrog impl", {
  parameters <- read_parameters(test_path("testdata/adult_parms_full.h5"))
  actual <- run_model(parameters)

  expected <- read_parameters(test_path("testdata/leapfrog_fit_full.h5"))

  expect_equal(actual$p_totpop, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(actual$births, expected$births)
  expect_equal(actual$p_background_deaths_totpop, expected$natdeaths)
  expect_equal(actual$p_hivpop, expected$hivpop1)
  expect_equal(actual$p_background_deaths_hivpop, expected$natdeaths_hivpop)
  expect_equal(actual$h_hivpop, expected$hivstrat_adult)
  expect_equal(actual$h_artpop, expected$artstrat_adult)
  expect_equal(actual$h_hiv_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$p_infections, expected$infections)
  expect_equal(actual$h_hiv_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$h_art_initiation, expected$artinit)
  expect_equal(actual$p_hiv_deaths, expected$hivdeaths)
})

test_that("model agrees with leapfrog impl", {
  parameters <- read_parameters(test_path("testdata/adult_parms_coarse.h5"))
  actual <- run_model(parameters, "HivCoarseAgeStratification")

  expected <- read_parameters(test_path("testdata/leapfrog_fit_coarse.h5"))

  expect_equal(actual$p_totpop, expected$totpop1)
  ## expected births doesn't have dim attribute so drop it for tests
  expect_equal(actual$births, expected$births)
  expect_equal(actual$p_background_deaths_totpop, expected$natdeaths)
  expect_equal(actual$p_hivpop, expected$hivpop1)
  expect_equal(actual$p_background_deaths_hivpop, expected$natdeaths_hivpop)
  expect_equal(actual$h_hivpop, expected$hivstrat_adult)
  expect_equal(actual$h_artpop, expected$artstrat_adult)
  expect_equal(actual$h_hiv_deaths_no_art, expected$aidsdeaths_noart)
  expect_equal(actual$p_infections, expected$infections)
  expect_equal(actual$h_hiv_deaths_art, expected$aidsdeaths_art)
  expect_equal(actual$h_art_initiation, expected$artinit)
  expect_equal(actual$p_hiv_deaths, expected$hivdeaths)
})
