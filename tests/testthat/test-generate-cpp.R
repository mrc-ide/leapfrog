test_that("can generate output parsing", {
  t <- tempfile()
  generate_build_r_output(t)

  result <- readLines(t)
  expect_true(any(grepl("// ========= DO NOT EDIT =========", result)))
  expect_true(any(grepl(paste0(
    "Rcpp::NumericVector r_total_population",
    "\\(ss\\.age_groups_pop \\* ss\\.num_genders \\* output_years\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "r_total_population\\.attr\\(\\\"dim\\\"\\) = Rcpp::NumericVector::",
    "create\\(ss\\.age_groups_pop, ss\\.num_genders, output_years\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "std::copy_n\\(state\\.total_population\\.data\\(\\), ",
    "state.total_population.size\\(\\), REAL\\(r_total_population\\)\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "Rcpp::_\\[\\\"total_population\\\"\\] = r_total_population,"),
    result)))
})
