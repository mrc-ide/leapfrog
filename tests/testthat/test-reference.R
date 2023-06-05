test_that("model agrees with leapfrog impl", {
  demp <- readRDS(test_path("testdata/demographic_projection_object.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters.rds"))
  actual <- run_base_model(demp, parameters, NULL, NULL)

  expected <- readRDS(test_path("testdata/leapfrog_fit.rds"))

  ## We're only reporting out last year atm so check final time point agrees
  expect_equal(actual$total_population, expected$totpop1[, , 61])
  expect_equal(actual$births, expected$births[61])
  expect_equal(actual$natural_deaths, expected$natdeaths[, , 61])
  expect_equal(actual$hiv_population, expected$hivpop1[, , 61])
  expect_equal(actual$hiv_natural_deaths, expected$natdeaths_hivpop[, , 61])
  expect_equal(actual$hiv_strat_adult, expected$hivstrat_adult[, , , 61])
  expect_equal(actual$art_strat_adult, expected$artstrat_adult[, , , 61])
  expect_equal(actual$aids_deaths_no_art, expected$aidsdeaths_noart[, , , 61])
  expect_equal(actual$infections, expected$infections[, , 61])
  expect_equal(actual$aids_deaths_art, expected$aidsdeaths_art[, , , , 61])
  expect_equal(actual$art_initiation, expected$artinit[, , , 61])
  expect_equal(actual$hiv_deaths, expected$hivdeaths[, , 61])
})
