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
     ## using as.integer() removes any dimensions. *1 to preserve shape and convert to int
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

assert_enum <- function(x, values, name = deparse(substitute(x))) {
  if (!(x %in% values)) {
    value_text <- paste(sprintf("'%s'", values), collapse = ", ")
    stop(sprintf("%s must be one of %s, got '%s'", name, value_text, x),
         call. = FALSE)
  }
  invisible(TRUE)
}

vcapply <- function(X, FUN, ...) { # nolint
  vapply(X, FUN, character(1), ...)
}
