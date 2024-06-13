#!/usr/bin/env Rscript

message("* Generating R interfaces")
## Source files manually so that if we break our generated code
## we can still regenerate the interfaces without compiling
root <- here::here()
source(file.path(root, "R/util.R"))
source(file.path(root, "R/util_assert.R"))
source(file.path(root, "R/generate_cpp_validation.R"))
source(file.path(root, "R/generate_cpp.R"))

frogger_file <- function(..., mustWork = TRUE) {
  path <- file.path(root, "inst", ...)
  if (mustWork && !file.exists(path)) {
    stop(sprintf("File at path %s does not exist", path))
  }
  path
}

generate_input_interface(dest = frogger_file("include/generated/model_input.hpp"))
generate_output_interface(dest = frogger_file("include/generated/model_output.hpp"))
generate_parameter_types(dest = frogger_file("include/generated/parameter_types.hpp"))
generate_state_types(dest = frogger_file("include/generated/state_types.hpp"))
generate_state_saver_types(dest = frogger_file("include/generated/state_saver_types.hpp"))
