#' Generate C++ for passing output from model fit into R
#'
#' This generates using metadata from `inst/cpp_generation/model_output.csv`
#'
#' @param dest The destination to write generated code to.
#'
#' @return Nothing, called to generate code in src dir
#' @keywords internal
generate_output_interface <- function(dest) {
  template_path <- frogger_file("cpp_generation/model_output.hpp.in")
  template <- readLines(template_path)
  output_file <- "model_output.csv"
  outputs <- utils::read.csv(frogger_file("cpp_generation", output_file),
                             colClasses = "character")

  validate_dimensions_columns(colnames(outputs), output_file)

  parsed_outputs <- lapply(seq_len(nrow(outputs)), function(row_num) {
    row <- outputs[row_num, ]
    ## When reading csv in excel the header column is included in count
    csv_row_num <- row_num + 1
    validate_and_parse_output(as.list(row), output_file, csv_row_num)
  })

  output_sections <- group_list_of_lists(parsed_outputs, "output_when")

  default_section <- output_sections[names(output_sections) == ""][[1]]
  ## TODO: don't hardcode the struct section here, remove this after
  ## we are generating out struct types and we know which output
  ## belongs to which struct
  default_data <- generate_output_section(default_section, "base")
  default_data <- paste_lines(default_data)
  build_list <- generate_build_list(default_section)
  build_list <- paste_lines(build_list)

  remaining_sections <- output_sections[names(output_sections) != ""]
  optional_data <- Map(generate_optional_output_section,
                       remaining_sections, names(remaining_sections))
  optional_data <- paste_lines(unlist(optional_data))

  header <- generate_header(basename(template_path))

  generated_code <- generate_cpp(template)
  writeLines(generated_code, dest)
  invisible(dest)
}

generate_output_section <- function(section, struct_name) {
  unpack <- generate_unpack_state_space(section)
  initialise_memory <- generate_initialise_r_memory(section)
  set_r_dimensions <- generate_set_r_dimensions(section)
  copy_data <- generate_copy_data(section, struct_name)
  c(unpack, initialise_memory, set_r_dimensions, copy_data)
}

generate_optional_output_section <- function(section, output_when) {
  condition <- sprintf("  if constexpr (%s) {", output_when)
  data <- generate_output_section(section, "children")
  push_to_list <- generate_push_to_list(section)
  c(condition, paste0("  ", c(data, push_to_list)), "  }")
}

generate_unpack_state_space <- function(outputs) {
  used <- unique(unlist(lapply(outputs, "[[", "parsed_dims")))
  ss <- used[grepl("^\\w+\\.\\w+$", used)]
  state_space_groups <- unique(vcapply(strsplit(ss, "\\."), "[[", 1))
  non_ss_names <- c("ss", "data", "options")
  state_space_groups <- state_space_groups[
    !(state_space_groups %in% non_ss_names)]
  vcapply(state_space_groups, function(group) {
    sprintf("  constexpr auto %s = ss.%s;", group, group)
  })
}

generate_initialise_r_memory <- function(outputs) {
  vcapply(outputs, function(output) {
    dimensions <- paste(output$parsed_dims, collapse = " * ")
    sprintf("  Rcpp::NumericVector r_%s(%s);", output$r_name, dimensions)
  })
}

generate_set_r_dimensions <- function(outputs) {
  vcapply(outputs, function(output) {
    dimensions <- paste(output$parsed_dims, collapse = ", ")
    sprintf("  r_%s.attr(\"dim\") = Rcpp::NumericVector::create(%s);",
            output$r_name, dimensions)
  })
}

generate_copy_data <- function(outputs, struct_name) {
  vcapply(outputs, function(output) {
    sprintf("  std::copy_n(state.%s.%s.data(), state.%s.%s.size(), %s(r_%s));",
            struct_name, output$cpp_name, struct_name, output$cpp_name,
            output$r_type, output$r_name)
  })
}

generate_build_list <- function(outputs) {
  initialise_list <- c(
    sprintf("  Rcpp::List ret(%s);", length(outputs)),
    sprintf("  Rcpp::CharacterVector names(%s);", length(outputs)))
  items <- vcapply(seq_along(outputs), function(i) {
    sprintf("  names[%s] = \"%s\";\n  ret[%s] = r_%s;",
            i - 1, outputs[[i]]$r_name, i - 1, outputs[[i]]$r_name)
  })
  set_names <- "  ret.attr(\"names\") = names;"
  c(initialise_list, items, set_names)
}

generate_push_to_list <- function(outputs) {
  vcapply(outputs, function(output) {
    sprintf("  ret.push_back(r_%s, \"%s\");", output$r_name, output$r_name)
  })
}

#' Generate C++ for passing input data in model fit
#'
#' This generates using metadata from `inst/cpp_generation/model_input.csv`
#'
#' @param dest The destination to write generated code to.
#' @param input_csv Path to the csv of model inputs.
#'
#' @return Nothing, called to generate code in src dir
#' @keywords internal
generate_input_interface <- function(
    dest, input_csv = frogger_file("cpp_generation/model_input.csv")) {

  template_path <- frogger_file("cpp_generation/model_input.hpp.in")
  template <- readLines(template_path)
  input_file <- basename(input_csv)
  inputs <- utils::read.csv(input_csv, colClasses = "character")

  validate_dimensions_columns(colnames(inputs), input_file)

  parsed_inputs <- lapply(seq_len(nrow(inputs)), function(row_num) {
    row <- inputs[row_num, ]
    ## When reading csv in excel the header column is included in count
    csv_row_num <- row_num + 1
    validate_and_parse_input(as.list(row), input_file, csv_row_num)
  })

  input_sections <- group_list_of_lists(parsed_inputs, "input_when")

  default_section <- input_sections[names(input_sections) == ""][[1]]
  default_data <- generate_input_section(default_section)
  default_data <- paste_lines(default_data)

  structs <- generate_struct_instantiation(default_section)
  structs <- paste_lines(structs)

  remaining_sections <- input_sections[names(input_sections) != ""]
  optional_data <- Map(generate_optional_input_section,
                       remaining_sections, names(remaining_sections))
  optional_data <- paste_lines(unlist(optional_data))

  header <- generate_header(basename(template_path))

  generated_code <- generate_cpp(template)
  writeLines(generated_code, dest)
  invisible(dest)
}


generate_input_section <- function(section) {
  from_r <- vlapply(section, function(input) {
    is_set(input$r_name)
  })
  from_value <- section[!from_r]
  from_r <- section[from_r]

  unpack <- generate_unpack_state_space(section)
  r_parse_data <- generate_input_from_r(from_r)
  value_data <- generate_input_from_value(from_value)
  c(unpack, r_parse_data, value_data)
}

generate_optional_input_section <- function(section, input_when) {
  condition <- sprintf("  if constexpr (%s) {", input_when)
  data <- generate_input_section(section)
  structs <- generate_struct_instantiation(section)
  return_statement <- generate_return()
  c(condition, paste0("  ", c(data, structs, return_statement)), "  }")
}

generate_struct_instantiation <- function(inputs) {
  inputs_by_struct <- get_data_by_struct(inputs)
  struct_text <- lapply(inputs_by_struct, generate_struct)
  unlist(struct_text)
}

generate_struct <- function(struct_inputs) {
  c(
    sprintf("  const leapfrog::%s<real_type> %s_params = {",
            struct_inputs[[1]]$struct,
            camel_to_snake(struct_inputs[[1]]$struct)),
    sprintf("      %s,", vcapply(struct_inputs, "[[", "cpp_name")),
    "  };"
  )
}

generate_input_from_r <- function(inputs) {
  vcapply(inputs, function(input) {
    dimensions <- paste(input$parsed_dims, collapse = ", ")
    if (input$dims == 1 && dimensions[1] == 1) {
      return(generate_length1_input(input))
    }
    rhs <- sprintf("parse_data<%s>(data, \"%s\", %s)",
                   input$type, input$r_name, dimensions)
    tensor_type <- "leapfrog::TensorMap"
    if (!is.null(input$convert_0_based) && input$convert_0_based) {
      rhs <- sprintf("convert_0_based<%s>(%s)", input$dims, rhs)
      ## Must be a tensor otherwise the convert_0_based will
      ## modify the underlying R data which we do not want to do
      tensor_type <- "leapfrog::Tensor"
    }
    lhs <- sprintf("  const %s%s<%s> %s",
                   tensor_type, input$dims, input$type, input$cpp_name)
    paste0(lhs, " = ", rhs, ";")
  })
}

generate_length1_input <- function(input) {
  sprintf("  const %s %s = Rcpp::as<%s>(data[\"%s\"]);",
          input$type, input$cpp_name, input$type, input$r_name)
}

generate_input_from_value <- function(inputs) {
  vcapply(inputs, function(input) {
    dimensions <- paste(input$parsed_dims, collapse = ", ")
    declaration <- sprintf("  leapfrog::Tensor%s<%s> %s(%s);",
                           input$dims, input$type, input$cpp_name,
                           input$parsed_dims)
    set_value <- sprintf("  %s.setConstant(%s);",
                         input$cpp_name, input$value)
    paste0(declaration, "\n", set_value)
  })
}

generate_header <- function(source_file) {
  paste0(
    "// Generated by frogger: do not edit by hand\n",
    "// This file is automatically generated. Do not edit this file. If you ",
    sprintf("want to make changes\n// edit `%s` and run ", source_file),
    "`./scripts/generate` to regenerate.")
}

generate_return <- function() {
  ## TODO: Generate this when we have type definitions being generated
  c("  const leapfrog::ChildModelParameters<ModelVariant, real_type> child_model_params = {",
    "      children_params",
    "  };",
    "  return leapfrog::Parameters<ModelVariant, real_type> {",
    "      base_model_params,",
    "      child_model_params",
    "  };",
    "} else {",
    "  return leapfrog::Parameters<ModelVariant, real_type> {",
    "      base_model_params",
    "  };")
}

#' Generate C++ for input parameter types
#'
#' This generates using metadata from `inst/cpp_generation/model_input.csv`
#'
#' @param dest The destination to write generated code to.
#'
#' @return Nothing, called to generate code in src dir
#' @keywords internal
generate_parameter_types <- function(dest) {
  template_path <- frogger_file("cpp_generation/parameter_types.hpp.in")
  template <- readLines(template_path)
  input_file <- "model_input.csv"
  inputs <- utils::read.csv(frogger_file("cpp_generation", input_file),
                            colClasses = "character")

  validate_dimensions_columns(colnames(inputs), input_file)

  parsed_inputs <- lapply(seq_len(nrow(inputs)), function(row_num) {
    row <- inputs[row_num, ]
    ## When reading csv in excel the header column is included in count
    csv_row_num <- row_num + 1
    validate_and_parse_input(as.list(row), input_file, csv_row_num)
  })

  inputs_by_struct <- get_data_by_struct(parsed_inputs)

  struct_defs <- vcapply(inputs_by_struct, generate_struct_def)
  struct_defs <- paste_lines(struct_defs)

  header <- generate_header(basename(template_path))

  generated_code <- generate_cpp(template)
  writeLines(generated_code, dest)
  invisible(dest)
}

#' Organise the data into separate structs
#'
#' This takes a list of inputs or outputs and splits them by struct returning
#' the result as a named list of lists where names are the struct name and
#' list is the inputs/outputs which belong on that struct
#'
#' @param data Data related to model inputs or outputs
#'
#' @return
#' @keywords internal
get_data_by_struct <- function(data) {
  struct <- vcapply(data, "[[", "struct")
  structs <- unique(vcapply(data, "[[", "struct"))
  data_by_struct <- lapply(structs, function(struct_name) {
    data[struct == struct_name]
  })
  names(data_by_struct) <- structs
  data_by_struct
}

generate_struct_def <- function(inputs) {
  input_text <- vcapply(inputs, function(input) {
    if (input$dims == 1 && input$dim1 == 1) {
      type <- input$type
    } else if (input$convert_0_based || nzchar(input$value)) {
      type <- sprintf("Tensor%s<%s>", input$dims, input$type)
    } else {
      type <- sprintf("TensorMap%s<%s>", input$dims, input$type)
    }
    sprintf("  %s %s;", type, input$cpp_name)
  })
  paste0(
    "template<typename real_type>\n",
    sprintf("struct %s {\n", inputs[[1]]$struct),
    paste_lines(input_text),
    "\n};\n"
  )
}

generate_struct_instance <- function(inputs) {
  paste0(
    sprintf("  const leapfrog::%s<real_type> %s_params = {\n",
            inputs[[1]]$struct,
            to_lower_camel(inputs[[1]]$struct)),
    paste_lines(sprintf("    %s", vcapply(inputs, "[[", "cpp_name"))),
    "  \n};\n"
  )
}

#' Generate C++ for state types
#'
#' @param dest The destination to write generated code to.
#' @param output_csv Path to csv file with output specification in it
#'
#' @return Nothing, called to generate code in src dir
#' @keywords internal
generate_state_types <- function(
    dest,
    output_csv = frogger_file("cpp_generation/model_output.csv")) {
  template_path <- frogger_file("cpp_generation/state_types.hpp.in")
  template <- readLines(template_path)
  output_file <- basename(output_csv)
  outputs <- utils::read.csv(output_csv, colClasses = "character")

  validate_dimensions_columns(colnames(outputs), output_file)

  parsed_outputs <- lapply(seq_len(nrow(outputs)), function(row_num) {
    row <- outputs[row_num, ]
    ## When reading csv in excel the header column is included in count
    csv_row_num <- row_num + 1
    output <- validate_and_parse_output(as.list(row), output_file, csv_row_num)
    output$parsed_dims <- parse_state_dims(output$parsed_dims)
    output
  })

  outputs_by_struct <- get_data_by_struct(parsed_outputs)

  state_defs <- vcapply(outputs_by_struct, generate_state_def)
  state_defs <- paste_lines(state_defs)

  header <- generate_header(basename(template_path))

  generated_code <- generate_cpp(template)
  writeLines(generated_code, dest)
  invisible(dest)
}

#' Generate C++ for state saver types
#'
#' @param dest The destination to write generated code to.
#' @param output_csv Path to csv file with output specification in it
#'
#' @return Nothing, called to generate code in src dir
#' @keywords internal
generate_state_saver_types <- function(
    dest,
    output_csv = frogger_file("cpp_generation/model_output.csv")) {
  template_path <- frogger_file("cpp_generation/state_saver_types.hpp.in")
  template <- readLines(template_path)
  output_file <- "model_output.csv"
  outputs <- utils::read.csv(output_csv, colClasses = "character")

  validate_dimensions_columns(colnames(outputs), output_file)

  parsed_outputs <- lapply(seq_len(nrow(outputs)), function(row_num) {
    row <- outputs[row_num, ]
    ## When reading csv in excel the header column is included in count
    csv_row_num <- row_num + 1
    output <- validate_and_parse_output(as.list(row), output_file, csv_row_num)
    output
  })

  outputs_by_struct <- get_data_by_struct(parsed_outputs)

  output_state_defs <- vcapply(outputs_by_struct, generate_output_state_defs)
  output_state_def <- paste_lines(output_state_defs)

  header <- generate_header(basename(template_path))

  generated_code <- generate_cpp(template)
  writeLines(generated_code, dest)
  invisible(dest)
}

parse_state_dims <- function(dims) {
  ## For output types the last dimension should always be output_years
  ## the state stores the state for a single year so to get the dimensions
  ## for generating state types we need to drop this last dimension
  ## We also don't need the qualifying name i.e. the `base` part of `base.hTS`
  ## We've already validated that the last dimension is output_years
  dims <- dims[-length(dims)]
  vcapply(dims, function(dim) {
    strsplit(dim, "\\.")[[1]][[2]]
  }, USE.NAMES = FALSE)
}

generate_state_def <- function(outputs) {
  struct_members <- vcapply(outputs, generate_state_struct_member)
  reset_text <- generate_state_reset_function(outputs)
  struct_name <- outputs[[1]]$struct
  model_variant <- outputs[[1]]$model_variant

  if (model_variant == "ModelVariant") {
    struct_def <- paste0(
      sprintf("template<typename %s, typename real_type>\n", model_variant),
      sprintf("struct %sState {\n", struct_name)
    )
  } else {
    struct_def <- paste0(
      "template<typename real_type>\n",
      sprintf("struct %sState<%s, real_type> {\n", struct_name, model_variant)
    )
  }

  paste0(
    struct_def,
    paste_lines(struct_members),
    "\n\n",
    sprintf("  %sState(const Parameters<%s, real_type> &pars) {\n",
            struct_name, model_variant),
    "    reset();\n",
    "  }\n\n",
    paste_lines(reset_text),
    "\n};\n"
  )
}

generate_state_struct_member <- function(output) {
  if (output$dims == "1") {
    type <- "real_type"
  } else {
    sizes <- sprintf("%s<%s>", output$parsed_dims, output$model_variant)
    sizes <- paste(sizes, collapse = ", ")
    type <- sprintf("TensorFixedSize<real_type, Sizes<%s>>", sizes)
  }
  sprintf("  %s %s;", type, output$cpp_name)
}

generate_state_reset_function <- function(outputs) {
  reset_lines <- vcapply(outputs, function(output) {
    if (output$dims == 1) {
      sprintf("    %s = 0;", output$cpp_name)
    } else {
      sprintf("    %s.setZero();", output$cpp_name)
    }
  })
  c(
    "  void reset() {",
    reset_lines,
    "  }")
}

generate_output_state_defs <- function(outputs) {
  struct_members <- vcapply(outputs, generate_output_state_struct_member)
  initialiser_list <- vcapply(outputs, generate_output_state_initialiser_list)
  ctor_body <- vcapply(outputs, function(output) {
    sprintf("    %s.setZero();", output$cpp_name)
  })
  reset_text <- generate_state_reset_function(outputs)
  struct_name <- outputs[[1]]$struct
  model_variant <- outputs[[1]]$model_variant

  if (model_variant == "ModelVariant") {
    struct_def <- paste0(
      sprintf("template<typename %s, typename real_type>\n", model_variant),
      sprintf("struct %sOutputState {\n", struct_name)
    )
  } else {
    struct_def <- paste0(
      "template<typename real_type>\n",
      sprintf("struct %sOutputState<%s, real_type> {\n", struct_name, model_variant)
    )
  }

  paste0(
    struct_def,
    paste_lines(struct_members),
    "\n\n",
    sprintf("  %sOutputState(int output_years): \n", struct_name),
    paste(initialiser_list, collapse = ",\n"),
    " {\n",
    paste_lines(ctor_body),
    "\n  }\n",
    "};\n"
  )
}

generate_output_state_struct_member <- function(output) {
  type <- sprintf("Tensor%s<real_type>", output$dims)
  sprintf("  %s %s;", type, output$cpp_name)
}

generate_output_state_initialiser_list <- function(output) {
  dim_indentation <- "      "
  dims <- vcapply(output$parsed_dims, function(dim) {
    if (dim == "output_years") {
      paste0(dim_indentation, dim)
    } else {
      sprintf("%sStateSpace<%s>().%s",
              dim_indentation, output$model_variant, dim)
    }
  })
  dims <- paste(dims, collapse = ",\n")
  sprintf("    %s(\n%s\n    )", output$cpp_name, dims)
}

generate_cpp <- function(template) {
  glue::glue(paste_lines(template),
             .open = "{{", .close = "}}",
             .envir = parent.frame())
}

camel_to_snake <- function(x) {
  has_uppercase_letter <- grepl("[A-Z]", x)
  if (!has_uppercase_letter) {
    return(x)
  }
  output_string <- gsub("([A-Z])", "_\\1", x)
  ## If input string starts with a capital this will append _ at the front
  ## remove it in this case
  if (substr(output_string, 1, 1) == "_") {
    output_string <- substr(output_string, 2, nchar(output_string))
  }

  tolower(output_string)
}
