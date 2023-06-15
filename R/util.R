`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}

## Serialize an R array into a format that can be read as an Eigen::Tensor
## by deserialize_tensor
serialize_r_to_tensor <- function(data, path) {
  type <- "double"
  if (is.integer(data)) {
    type <- "int"
  }
  dims <- paste0(dim(data), collapse = ",")
  data <- paste0(data, collapse = ",")
  lines <- c(type, dims, data)
  writeLines(lines, path)
  path
}
