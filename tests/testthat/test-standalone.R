test_that("We can compile the standalone program", {
  skip_for_compilation()
  skip_on_os("windows")

  path_src <- frogger_file("fit_model")
  tmp <- tempfile()
  copy_directory(path_src, tmp)

  args <- dirname(frogger_file("include"))

  code <- withr::with_dir(
    tmp,
    system2("./configure", args, stdout = FALSE, stderr = FALSE))
  expect_equal(code, 0)
  code <- withr::with_dir(
    tmp,
    system2("make", stdout = FALSE, stderr = FALSE))
  expect_equal(code, 0)

  output <- tempfile()
  input <- frogger_file("fit_model/data")

  res <- system2(file.path(tmp, "fit_model"),
                 c(60, 10, input, output),
                 stdout = TRUE)
  expect_equal(res, c(paste0("Created output directory '", output, "'"),
                      "Fit complete"))

  expect_true(file.exists(file.path(output, "total_population")))
  expect_true(file.exists(file.path(output, "births")))
  expect_true(file.exists(file.path(output, "natural_deaths")))
  expect_true(file.exists(file.path(output, "hiv_population")))
  expect_true(file.exists(file.path(output, "hiv_natural_deaths")))
  expect_true(file.exists(file.path(output, "hiv_strat_adult")))
  expect_true(file.exists(file.path(output, "art_strat_adult")))
  expect_true(file.exists(file.path(output, "aids_deaths_no_art")))
  expect_true(file.exists(file.path(output, "infections")))
  expect_true(file.exists(file.path(output, "aids_deaths_art")))
  expect_true(file.exists(file.path(output, "art_initiation")))
  expect_true(file.exists(file.path(output, "hiv_deaths")))

  expected <- readRDS(test_path("testdata/leapfrog_fit.rds"))

  ## We're only reporting out last year atm so check final time point agrees
  ## We're also expecting some precision loss due to serialization so check up to some tolerance
  expect_equal(deserialize_tensor_to_r(file.path(output, "total_population")),
               expected$totpop1[, , 61],
               tolerance = 1e-5)
  expect_equal(as.numeric(readLines(file.path(output, "births"))),
               expected$births[61],
               tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "natural_deaths")),
                 expected$natdeaths[, , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "hiv_population")),
                 expected$hivpop1[, , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "hiv_natural_deaths")),
                 expected$natdeaths_hivpop[, , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "hiv_strat_adult")),
                 expected$hivstrat_adult[, , , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "art_strat_adult")),
                 expected$artstrat_adult[, , , , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "aids_deaths_no_art")),
                 expected$aidsdeaths_noart[, , , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "infections")),
                 expected$infections[, , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "aids_deaths_art")),
                 expected$aidsdeaths_art[, , , , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "art_initiation")),
                 expected$artinit[, , , 61],
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "hiv_deaths")),
                 expected$hivdeaths[, , 61],
                 tolerance = 1e-5)
})
