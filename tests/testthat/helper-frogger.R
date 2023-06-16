skip_for_compilation <- function() {
  testthat::skip_on_cran()
}

frogger_file <- function(path) {
  system.file(path, package = "frogger", mustWork = TRUE)
}

copy_directory <- function(src, as) {
  files <- dir(src, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  dir.create(as, FALSE, TRUE)
  ok <- file.copy(files, as, recursive = TRUE)
  if (!all(ok)) {
    stop("Error copying files")
  }
}
