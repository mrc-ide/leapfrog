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

  expect_true(file.exists(file.path(output, "hiv_population")))
  content <- readLines(file.path(output, "hiv_population"))
  expect_equal(content[1], "double")
  expect_equal(content[2], "81,2")
  hiv_pop <- array(unlist(read.csv(text = content[3], header = FALSE)),
                   unlist(read.csv(text = content[2], header = FALSE)))

  expected <- readRDS(test_path("testdata/leapfrog_fit.rds"))

  ## We're only reporting out last year atm so check final time point agrees
  ## We're also expecting some precision loss due to serialization so check up to some tolerance
  expect_equal(hiv_pop, expected$hivpop1[, , 61], tolerance = 1e-6, ignore_attr = TRUE)
})
