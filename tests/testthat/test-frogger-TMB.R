test_that("initial state set up works as expected", {
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  report <- test_tmb(demp, parameters, NULL, NULL, run_child_model = FALSE)
  print(report)
})
