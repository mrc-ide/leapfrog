## usethis namespace: start
#' @useDynLib frogger, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
TMB::compile("src/frogger_TMB.cpp")
NULL
