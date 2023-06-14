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

  res <- system2(file.path(tmp, "fit_model"),
                 stdout = TRUE)
  expect_equal(res, "Ran thing")
})
