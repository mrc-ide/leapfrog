.onUnload <- function(libpath) { # nolint # nocov start
  library.dynam.unload("frogger", libpath)
} # nocov end

## This mess is to avoid R CMD check NOTEs about "no visible binding for
## global variable". These NOTEs aren't that important as we're using
## lots of dplyr here, but they do obscure real issues in the CMD check
## output. So declare a list of global variables to avoid this.
utils::globalVariables(c("Value", "Sex", "Age", "Year"))

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
