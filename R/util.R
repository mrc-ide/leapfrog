`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}

#' Serialize an R array into a format that can be read as an Eigen::Tensor
#' by deserialize_tensor
#'
#' @param data The data to serialize, should be an array or vector
#' @param path File path to serialize data to
#'
#' @return Path to file
#'
#' @keywords internal
serialize_r_to_tensor <- function(data, path) {
  type <- "double"
  if (is.integer(data)) {
    type <- "int"
  } else if (is.logical(data)) {
    ## using as.integer() removes any dimensions.
    ## *1 to preserve shape and convert to int
    data <- data * 1
    type <- "int"
  }
  dims <- dim(data)
  if (is.null(dims)) {
    ## Assume 1 dimensional
    dims <- length(data)
  }
  dims <- paste0(dims, collapse = ",")
  data <- paste0(data, collapse = ",")
  lines <- c(type, dims, data)
  writeLines(lines, path)
  path
}

#' Deseralize an Eigen::Tensor into an R array
#'
#' @param path Path to serialized data
#'
#' @return An R array containing the deserialized data
#'
#' @keywords internal
deserialize_tensor_to_r <- function(path) {
  content <- readLines(path)
  array(as.numeric(strsplit(content[[3]], ",\\s*")[[1]]),
        as.numeric(strsplit(content[[2]], ",\\s*")[[1]]))
}

vcapply <- function(X, FUN, ...) { # nolint
  vapply(X, FUN, character(1), ...)
}

vlapply <- function(X, FUN, ...) { # nolint
  vapply(X, FUN, logical(1), ...)
}

frogger_file <- function(..., mustWork = TRUE) {
  system.file(..., package = "frogger", mustWork = mustWork)
}

is_unset <- function(x) {
  is.null(x) || is.na(x) || trimws(x) == ""
}

is_set <- function(x) {
  !is_unset(x)
}

format_vector <- function(vector) {
  paste(paste0("'", vector, "'"), collapse = ", ")
}


#' Take a list of lists and split into groups based on some property
#'
#' Every list must have this property otherwise this will fail
#'
#' @param list List of lists to split
#' @param on Name of property to split on
#'
#' @return List of lists of lists :O,
#' @keywords internal
#'
#' @examples
#' input_list <- list(
#'   list(id = 1, type = "foo"),
#'   list(id = 2, type = "foo"),
#'   list(id = 3, type = "bar"),
#'   list(id = 4, type = "bar")
#' )
#' group_list_of_lists(input_list, "type")
group_list_of_lists <- function(list, on) {
  property <- vcapply(list, "[[", on)
  types <- unique(property)
  out <- lapply(types, function(type) {
    list[property == type]
  })
  names(out) <- types
  out
}
