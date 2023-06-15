`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}

## Serialize an R array into a format that can be read as an Eigen::Tensor
## by deserialize_tensor
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
