## Script to compare delphi outputs to leapfrog ones

pjnz <- "C:/Users/Test/Downloads/demo_mwi2024_v6.43.PJNZ"
param_dirs <- c("C:/Users/Test/Downloads/demProjParams",
                "C:/Users/Test/Downloads/hivAdultParams",
                "C:/Users/Test/Downloads/hivChildParams")
output_dirs <- c("C:/Users/Test/Downloads/demProjState",
                 "C:/Users/Test/Downloads/hivAdultState",
                 "C:/Users/Test/Downloads/hivChildState")

repo_root <- gert::git_find()
configs_dir <- file.path(repo_root, "cpp_generation/modelSchemas/configs")

to_lower_camel <- function(x) {
  sp <- strsplit(tolower(x), "_", TRUE)
  vapply(sp, function(x) {
    out <- x[1]
    if (length(x) > 1) {
      upper <- vapply(x[-1], function(k) paste0(toupper(substr(k, 1, 1)), substr(k, 2, nchar(k))), character(1))
      out <- paste(c(out, upper), collapse = "")
    }
    out
  }, character(1))
}

to_snake_case <- function(x) {
  tolower(gsub("(?<=[a-z0-9])(?=[A-Z])", "_", x, perl = TRUE))
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
  setNames(all_params, to_lower_camel(names(all_params)))
}

name_mapping <- build_name_mapping()

pkgload::load_all()
source('./scripts/spectrum_inputs_paeds.R')
source('./scripts/read_spectrum.R')
demp <- prepare_leapfrog_demp(pjnz)
proj <- prepare_leapfrog_projp(pjnz)
proj <- prepare_hc_leapfrog_projp(pjnz, proj)
params <- c(demp, proj)
params$mat_prev_input <- rep(FALSE, 61)
params$basepop <- params$basepop[, , 1]
# scale_cd4_mort is always true from Delphi so force true here
params$scale_cd4_mort <- TRUE
expected <- process_parameters(params, "ChildModel")

delphi_params <- list.files(param_dirs, full.names = TRUE)
delphi_params <- setNames(delphi_params, basename(delphi_params))
actual <- lapply(delphi_params, deserialize_tensor_to_r)

## There are some inputs which you can supply from R to use direct input instead
## of using data from the adult model. For delphi we'll always want to use
## the adult model, so we are expecting differences in these values
dont_compare <- c("matPrevInput", "propLt200",
                  "propGte350", "matHivBirths",
                  "adultFemaleInfections", "totalBirths")

compare <- names(actual)[!(names(actual) %in% dont_compare)]

for (param in compare) {
  message("Comparing ", param)
  r_name <- name_mapping[[param]]
  compare_to <- expected[[r_name]]
  if (is.logical(compare_to)) {
    compare_to <- as.numeric(compare_to)
  }
  testthat::expect_equal(actual[[param]], compare_to,
                         ignore_attr = TRUE, label = param,
                         expected.label = r_name)
}

expected_result <- run_model(params, "ChildModel")
delphi_output <- list.files(output_dirs, full.names = TRUE)
delphi_output <- setNames(delphi_output, basename(delphi_output))
actual_result <- lapply(delphi_output, deserialize_tensor_to_r)

for (state in names(actual_result)) {
  message("Comparing ", state)
  r_name <- to_snake_case(state)
  compare_to <- expected_result[[r_name]]
  if (is.logical(compare_to)) {
    compare_to <- as.numeric(compare_to)
  }
  testthat::expect_equal(actual_result[[state]], compare_to,
                         ignore_attr = TRUE, label = state,
                         expected.label = r_name,
                         tolerance = 1e-7)
}
