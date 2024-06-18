testthat::skip("Standalone tests currently broken")

test_that("We can compile the standalone program", {
  skip_for_compilation()
  skip_on_os("windows")
  testthat::skip_if_not(
    file.exists(frogger_file("standalone_model/data/adult_data",
                             mustWork = FALSE)),
    paste("Skipping as test data not extracted, run ",
          "'./inst/standalone_model/extract_data' first."))

  path_src <- frogger_file("standalone_model")
  tmp <- tempfile()
  copy_directory(path_src, tmp)

  # Remove files if they have already been compiled to ensure we are always
  # compiling in tihs test to ensure linking to latest frogger code
  unlink(file.path(tmp, "Makefile"))
  unlink(file.path(tmp, "simulate_model"))

  args <- c(dirname(frogger_file("include")), find.package("RcppEigen"))

  code <- withr::with_dir(
    tmp,
    system2("./configure", args, stdout = FALSE, stderr = FALSE)
  )
  expect_equal(code, 0)
  code <- withr::with_dir(
    tmp,
    system2("make", "-B", stdout = FALSE, stderr = FALSE)
  )
  expect_equal(code, 0)

  output <- tempfile()
  input <- frogger_file("standalone_model/data/adult_data")

  res <- system2(file.path(tmp, "simulate_model"),
    c(61, 10, input, output),
    stdout = TRUE
  )
  expect_equal(res, c(
    paste0("Created output directory '", output, "'"),
    "Fit complete"
  ))

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

  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))
  expected <- run_model(demp, parameters, 1970:2030L, 10L,
                        run_child_model = FALSE)

  ## There is some precision loss in standalone data due to serialization
  ## so check up to appropriate tolerance
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_total_pop")),
    expected$p_total_pop,
    tolerance = 1e-5
  )
  expect_equal(as.numeric(deserialize_tensor_to_r(file.path(output, "births"))),
    expected$births,
    tolerance = 1e-5,
    ignore_attr = TRUE
  )
  expect_equal(
    deserialize_tensor_to_r(file.path(output, "p_total_pop_natural_deaths")),
    expected$p_total_pop_natural_deaths,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_hiv_pop")),
    expected$p_hiv_pop,
    tolerance = 1e-5
  )
  expect_equal(
    deserialize_tensor_to_r(file.path(output, "p_hiv_pop_natural_deaths")),
    expected$p_hiv_pop_natural_deaths,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_hiv_adult")),
    expected$h_hiv_adult,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_art_adult")),
    expected$h_art_adult,
    tolerance = 1e-5
  )
  expect_equal(
    deserialize_tensor_to_r(file.path(output, "h_hiv_deaths_no_art")),
    expected$h_hiv_deaths_no_art,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_infections")),
    expected$p_infections,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_hiv_deaths_art")),
    expected$h_hiv_deaths_art,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "h_art_initiation")),
    expected$h_art_initiation,
    tolerance = 1e-5
  )
  expect_equal(deserialize_tensor_to_r(file.path(output, "p_hiv_deaths")),
    expected$p_hiv_deaths,
    tolerance = 1e-5
  )
})
