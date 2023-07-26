test_that("can generate output parsing", {
  t <- tempfile()
  generate_output_interface(t)

  result <- readLines(t)
  expect_true(any(grepl("// ========= DO NOT EDIT =========", result)))
  expect_true(any(grepl(paste0(
    "Rcpp::NumericVector r_total_population",
    "\\(age_groups_pop \\* num_genders \\* output_years\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "r_total_population\\.attr\\(\\\"dim\\\"\\) = Rcpp::NumericVector::",
    "create\\(age_groups_pop, num_genders, output_years\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "std::copy_n\\(state\\.total_population\\.data\\(\\), ",
    "state.total_population.size\\(\\), REAL\\(r_total_population\\)\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "Rcpp::_\\[\\\"total_population\\\"\\] = r_total_population,"),
    result)))
})

test_that("can generate input parsing", {
  t <- tempfile()
  generate_input_interface(t)

  result <- readLines(t)
  expect_true(any(grepl("// ========= DO NOT EDIT =========", result)))
  expect_true(any(grepl(paste0(
    "const leapfrog::TensorMap2<real_type> base_pop = parse_data<real_type>",
    "\\(data, \"basepop\", age_groups_pop, num_genders\\);"),
    result)))
  expect_true(any(grepl(paste0(
    "const leapfrog::TensorMap1<int> artcd4elig_idx = convert_base<1>\\(",
    "parse_data<int>\\(data, \"artcd4elig_idx\", proj_years \\+ 1\\)\\);"),
    result)))
  expect_true(any(grepl(
    "leapfrog::Tensor1<real_type> h_art_stage_dur\\(treatment_stages - 1\\);",
    result)))
  expect_true(any(grepl(
    "h_art_stage_dur.setConstant\\(0.5\\);",
    result)))
})

test_that("generated files are up to date", {
  target_input_file <- frogger_file("r_interface/model_input.hpp")
  t_input <- tempfile()
  generate_input_interface(t_input)
  expect_identical(
    readLines(t_input), readLines(target_input_file),
    info = paste0("Your interface is out of date, regenerate by running ",
                  "./scripts/generate")
  )

  target_output_file <- frogger_file("r_interface/model_output.hpp")
  t_output <- tempfile()
  generate_output_interface(t_output)
  expect_identical(
    readLines(t_output), readLines(target_output_file),
    info = paste0("Your interface is out of date, regenerate by running ",
                  "./scripts/generate")
  )
})
