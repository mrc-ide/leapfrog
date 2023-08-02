test_that("We can compile the standalone program", {
  skip_for_compilation()
  skip_on_os("windows")

  path_src <- frogger_file("fit_model")
  tmp <- tempfile()
  copy_directory(path_src, tmp)

  args <- c(dirname(frogger_file("include")), find.package("RcppEigen"))

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
                 c(60, 10, "false", input, output),
                 stdout = TRUE)
  expect_equal(res, c(paste0("Created output directory '", output, "'"),
                      "Fit complete"))

  expect_true(file.exists(file.path(output, "p_total_pop")))
  expect_true(file.exists(file.path(output, "births")))
  expect_true(file.exists(file.path(output, "p_total_pop_natural_deaths")))
  expect_true(file.exists(file.path(output, "p_hiv_pop")))
  expect_true(file.exists(file.path(output, "p_hiv_pop_natural_deaths")))
  expect_true(file.exists(file.path(output, "h_hiv_adult")))
  expect_true(file.exists(file.path(output, "h_art_adult")))
  expect_true(file.exists(file.path(output, "h_hiv_deaths_no_art")))
  expect_true(file.exists(file.path(output, "p_infections")))
  expect_true(file.exists(file.path(output, "h_hiv_deaths_art")))
  expect_true(file.exists(file.path(output, "h_art_initiation")))
  expect_true(file.exists(file.path(output, "p_hiv_deaths")))

  expected <- readRDS(test_path("testdata/leapfrog_fit.rds"))

  ## We're only reporting out last year atm so check final time point agrees
  ## We're also expecting some precision loss due to serialization so check up to some tolerance
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_total_pop")),
               expected$totpop1,
               tolerance = 1e-5)
  expect_equal(as.numeric(deserialize_tensor_to_r(file.path(output, "births"))),
               expected$births,
               tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_total_pop_natural_deaths")),
                 expected$natdeaths,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_hiv_pop")),
                 expected$hivpop1,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_hiv_pop_natural_deaths")),
                 expected$natdeaths_hivpop,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_hiv_adult")),
                 expected$hivstrat_adult,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_art_adult")),
                 expected$artstrat_adult,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_hiv_deaths_no_art")),
                 expected$aidsdeaths_noart,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_infections")),
                 expected$infections,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_hiv_deaths_art")),
                 expected$aidsdeaths_art,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_art_initiation")),
                 expected$artinit,
                 tolerance = 1e-5)
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_hiv_deaths")),
                 expected$hivdeaths,
                 tolerance = 1e-5)
})
