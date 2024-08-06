#' Generate C++ for passing input data in model fit
#'
#' This generates using metadata from `inst/cpp_generation/model_input.csv`
#'
#' @param dest The destination to write generated code to.
#' @param input_csv Path to the csv of model inputs.
#'
#' @return Nothing, called to generate code in src dir
#' @keywords internal
generate_input_interface_tmb <- function(
    dest, input_csv = frogger_file("cpp_generation/model_input.csv")) {
  template_path <- frogger_file("cpp_generation/model_input_TMB_new.hpp.in")
  template <- readLines(template_path)
  input_sections <- read_input_csv(input_csv)
  input_sections <- lapply(input_sections, add_member_types)

  tmb_data_members <- generate_tmb_members(input_sections)
  constructors <- generate_constructors(input_sections)

  section_names <- names(input_sections)
  full_input_sections <- lapply(section_names, function(x) get_relevant_input_sections(x, input_sections))
  names(full_input_sections) <- section_names

  if_bodies <- list()
  else_body <- list()
  for (name in section_names) {
    if (name == "") {
      section <- full_input_sections[[1]]
    } else {
      section <- full_input_sections[[name]]
    }
    section <- add_member_types(section)
    eigen_tensor_setup <- generate_eigen_tensor_setup(section)

    structs <- generate_struct_instantiation_tmb(section)

    section_return <- generate_return_tmb(name != "")

    section_body <- unlist(c(eigen_tensor_setup, structs, section_return))
    if (name == "") {
      else_body <- unlist(add_indentation(section_body, 1))
    } else {
      if_body <- c(sprintf("  if constexpr (%s) {", name), "    constexpr auto children = ss.children;", unlist(add_indentation(section_body, 2)), "  }")
      if_bodies <- c(if_bodies, if_body)
    }
  }

  function_body <- paste_lines(c(if_bodies, else_body))

  # eigen_tensor_setup <- generate_eigen_tensor_setup(input_sections)

  # structs <- generate_struct_instantiation_tmb(input_sections)
  # structs <- paste_lines(structs)

  # return_struct <- generate_return_tmb()
  # return_struct <- paste_lines(return_struct)

  generate_cpp(template, dest, basename(template_path))
  invisible(dest)
}

generate_tmb_members <- function(input_sections) {
  members <- list()
  for (section in input_sections) {
    section_members <- lapply(section, get_tmb_type_for_param)
    section_members <- Filter(Negate(is.null), section_members)

    members <- c(members, add_indentation(section_members, 1))
  }
  paste_lines(members)
}

get_tmb_type_for_param <- function(param) {
  structure_type <- param$member_type$structure_type
  data_type <- param$member_type$data_type
  cpp_member_type <- ""
  if (structure_type == "number") {
    cpp_member_type <- sprintf("%s", data_type)
  } else if (structure_type == "convert_0_based_tensor" || structure_type == "tensor_map") {
    cpp_member_type <- sprintf("array<%s>", data_type)
  } else {
    return(NULL)
  }
  sprintf("%s %s;", cpp_member_type, param$cpp_name)
}

generate_constructors <- function(input_sections) {
  constructors <- list()
  for (section in input_sections) {
    section_constructors <- lapply(section, get_section_constructor)
    if (section[[1]]$input_when != "") {
      section_constructors <- c(sprintf("if constexpr (%s) {", section[[1]]$input_when), add_indentation(section_constructors, 1), "}")
    }
    constructors <- c(constructors, section_constructors)
  }

  paste_lines(add_indentation(unlist(constructors), 2))
}

get_section_constructor <- function(param) {
  structure_type <- param$member_type$structure_type
  data_type <- param$member_type$data_type
  if (structure_type == "number") {
    if (data_type == "int") {
      sprintf("%s = CppAD::Integer(tmbutils::asVector<Type>(REAL(%s), 1)[0]);", param$cpp_name, data_from_sexp(param$r_name))
    } else {
      sprintf("%s = tmbutils::asVector<Type>(REAL(%s), 1)[0];", param$cpp_name, data_from_sexp(param$r_name))
    }
  } else if (structure_type == "convert_0_based_tensor") {
    sprintf("%s = tmbutils::asArray<%s>(%s);", param$cpp_name, data_type, data_from_sexp(param$r_name))
  } else if (structure_type == "tensor_map") {
    sprintf("%s = tmbutils::asArray<%s>(%s);", param$cpp_name, data_type, data_from_sexp(param$r_name))
  } else {
    NULL
  }
}

generate_eigen_tensor_setup <- function(section) {
  section_tmb_param <- lapply(section, get_section_tmb_param)
  unlist(section_tmb_param)
}

get_section_tmb_param <- function(param) {
  structure_type <- param$member_type$structure_type
  data_type <- param$member_type$data_type
  dimensions <- paste(static_cast_to_int(param$parsed_dims), collapse = ", ")
  if (structure_type == "number") {
    sprintf("const %s %s = %s;", data_type, param$cpp_name, get_tmb_data_elem(param$cpp_name))
  } else if (structure_type == "convert_0_based_tensor") {
    e_type <- sprintf("Eigen::Tensor<%s, %s>", data_type, param$dims)
    sprintf("const %s %s = convert_0_based_tmb_new<%s>(%s, %s);", e_type, param$cpp_name, data_type, get_tmb_data_elem(param$cpp_name), dimensions)
  } else if (structure_type == "tensor_map") {
    e_type <- sprintf("Eigen::TensorMap<const Eigen::Tensor<%s, %s>>", data_type, param$dims)
    sprintf("const %s %s = %s(%s.data(), %s);", e_type, param$cpp_name, e_type, get_tmb_data_elem(param$cpp_name), dimensions)
  } else if (structure_type == "constant_value_tensor") {
    e_type <- sprintf("Eigen::Tensor<%s, %s>", data_type, param$dims)
    list(
      sprintf("%s %s(%s);", e_type, param$cpp_name, dimensions),
      sprintf("%s.setConstant(%s);", param$cpp_name, param$value)
    )
  } else {
    NULL
  }
}

static_cast_to_int <- function(parsed_dims) {
  lapply(parsed_dims, function(x) {
    sprintf("static_cast<int>(%s)", x)
  })
}

get_tmb_data_elem <- function(name) {
  sprintf("tmb_data.%s", name)
}

get_member_type <- function(type, structure_type) {
  list(data_type = ifelse(type == "real_type", "Type", type), structure_type = structure_type)
}

add_member_types <- function(section) {
  lapply(section, function(x) {
    if (!is_set(x$r_name)) {
      structure_type <- "constant_value_tensor"
    } else if (x$dims == 1 && x$dim1 == 1) {
      structure_type <- "number"
    } else if (!is.null(x$convert_0_based) && x$convert_0_based) {
      structure_type <- "convert_0_based_tensor"
    } else {
      structure_type <- "tensor_map"
    }
    x$member_type <- get_member_type(x$type, structure_type)
    x
  })
}

add_indentation <- function(string_list, n_indents) {
  lapply(string_list, function(s) {
    paste0(strrep("  ", n_indents), s)
  })
}

generate_struct_instantiation_tmb <- function(section) {
  inputs_by_struct <- group_list_of_lists(section, "struct")
  struct_text <- lapply(inputs_by_struct, generate_struct_tmb)
  unlist(struct_text)
}

generate_struct_tmb <- function(struct_inputs) {
  c(
    sprintf(
      "const leapfrog::%s<Type> %s_params = {",
      struct_inputs[[1]]$struct,
      camel_to_snake(struct_inputs[[1]]$struct)
    ),
    add_indentation(sprintf("%s,", vcapply(struct_inputs, "[[", "cpp_name")), 1),
    "};"
  )
}

data_from_sexp <- function(name) {
  sprintf("getListElement(%s, \"%s\")", "data", name)
}

generate_return_tmb <- function(is_child_model) {
  ## TODO: Generate this when we have type definitions being generated
  base_model_params <- c(
    "const leapfrog::BaseModelParameters<Type> base_model_params = {",
    "  options,",
    "  demography_params,",
    "  incidence_params,",
    "  natural_history_params,",
    "  art_params",
    "};"
  )
  if (is_child_model) {
    c(
      base_model_params,
      "const leapfrog::ChildModelParameters<ModelVariant, Type> child_model_params = {",
      "  children_params",
      "};",
      "return leapfrog::Parameters<ModelVariant, Type> {",
      "  base_model_params,",
      "  child_model_params",
      "};"
    )
  } else {
    c(
      base_model_params,
      "return leapfrog::Parameters<ModelVariant, Type> {",
      "  base_model_params",
      "};"
    )
  }
}

get_relevant_input_sections <- function(input_when, input_sections) {
  if (input_when == "") {
    input_sections[[1]]
  } else {
    c(input_sections[[1]], input_sections[[2]])
  }
}
