assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
}

assert_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) {
    stop(sprintf("'%s' must be character", name), call. = FALSE)
  }
}

assert_scalar_character <- function(x, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_character(x, name)
}

assert_enum <- function(x, values, name = deparse(substitute(x))) {
  if (!(x %in% values)) {
    value_text <- format_vector(values)
    stop(sprintf("%s must be one of %s, got '%s'", name, value_text, x),
         call. = FALSE)
  }
  invisible(TRUE)
}

assert_set <- function(x, name = deparse(substitute(x))) {
  assert_scalar_character(x)
  if (is_unset(x)) {
    stop(sprintf("%s must be a non-null, non-na, non-empty string", name))
  }
}

assert_one_is_set <- function(items, names, name = deparse(substitute(items))) {
  has_a_value <- !vlapply(items[names], is_unset)
  if (sum(has_a_value) > 1) {
    non_null <- format_vector(names(has_a_value))
    stop(sprintf("Items %s from '%s' all have values, only one should be set",
                 non_null, name))
  }
}
