#!/usr/bin/env Rscript

## Run this if you want to recreate the test data for the standalone model.
## You might do this is say we add some new required data or if the data we are testing with
## changes for any reason

repo_root <- gert::git_find()
configs_dir <- file.path(repo_root, "cpp_generation/modelSchemas/configs")
save_root_dir <- file.path(repo_root, "r-package/inst/standalone_model/data")

write_standalone_data <- function(input_data, data_map, dest) {
  t <- tempfile()
  dir.create(t)

  ## Basepop handled separately as we only save a slice of this out
  frogger:::serialize_r_to_tensor(input_data$basepop[, , 1],
                                  file.path(t, "base_pop"))
  for (r_name in setdiff(names(data_map), "basepop")) {
    if (!is.null(input_data[[r_name]])) {
      frogger:::serialize_r_to_tensor(input_data[[r_name]], file.path(t, data_map[[r_name]]))
    }
  }
  o <- zip::zip(dest,
                list.files(t, full.names = TRUE),
                include_directories = FALSE,
                root = save_root_dir,
                mode = "cherry-pick")
  file.path(save_root_dir, o)
}

build_name_mapping <- function() {
  cfgs <- list.files(configs_dir, full.names = TRUE)
  all_params <- lapply(cfgs, function(cfg) {
    params <- jsonlite::read_json(cfg)$pars$default
    lapply(params, function(param) {
      param$alias$r
    })
  })
  all_params <- unlist(all_params)
  stats::setNames(names(all_params), all_params)
}

prepare_input_data <- function(input_data) {
  # Convert indices to 0 based
  if ("artcd4elig_idx" %in% names(input_data)) {
    input_data[["artcd4elig_idx"]] <- input_data[["artcd4elig_idx"]] - 1L
  }
  if ("paed_art_elig_cd4" %in% names(input_data)) {
    input_data[["paed_art_elig_cd4"]] <- input_data[["paed_art_elig_cd4"]] - 1L
  }

  hTS <- dim(input_data[["artmx_timerr"]])[[1]]
  input_data$h_art_stage_dur <- rep(0.5, hTS - 1)
  input_data
}

## Map from name in parameters to name expected in C++ code
name_mapping <- build_name_mapping()

parameters <- frogger::read_hdf5_file(testthat::test_path("testdata/adult_parms.h5"))
input_data <- prepare_input_data(parameters)

path <- write_standalone_data(input_data, name_mapping, "adult_data.zip")
message(sprintf("Wrote adult test data to %s", path))

source(file.path("R", "spectrum_inputs.R"))
child <- frogger::read_hdf5_file(testthat::test_path("testdata/child_parms.h5"))
input_data <- prepare_input_data(child$parameters)

path <- write_standalone_data(input_data, name_mapping, "child_data.zip")
message(sprintf("Wrote child test data to %s", path))
