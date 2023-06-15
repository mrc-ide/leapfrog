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
  array(unlist(read.csv(text = content[3], header = FALSE), use.names = FALSE),
        unlist(read.csv(text = content[2], header = FALSE), use.names = FALSE))
}
