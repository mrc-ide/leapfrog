required_files <- c(
  "adult_parms_full.h5",
  "adult_parms_coarse.h5",
  "child_test_utils.rds",
  "child_parms_full.h5",
  "child_parms_coarse.h5"
)

required_file_paths <- unlist(lapply(required_files, function(f) test_path(sprintf("testdata/%s", f))))

if (!all(file.exists(required_file_paths))) {
  withr::local_options(
    list(error = NULL),
    .local_envir = teardown_env()
  )
  stop("
Ribbit?
       _   _
      (.)_(.)
   _ (   _   ) _
  / \\/`-----'\\/ \\
__\\ ( (     ) ) /__
)   /\\ \\._./ /\\   (
 )_/ /|\\   /|\\ \\_(

Oops, looks like you don't have the test data generated, please run scripts/create_test_data.R from the root of the package!

")
}
