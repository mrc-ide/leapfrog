#!/usr/bin/env Rscript
root_dir <- here::here()
parameters <- frogger::read_parameters(testthat::test_path("testdata/adult_parms.h5"))
out <- frogger::run_model(parameters)

output_dir <- file.path(root_dir, "..", "model-outputs")
if (!fs::dir_exists(output_dir)) {
  fs::dir_create(output_dir)
}

frogger:::save_hdf5_file(out, file.path(output_dir, "r-output.h5"))
