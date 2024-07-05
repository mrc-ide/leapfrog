test_that("initial state set up works as expected", {
  TMB::compile("src/frogger_TMB.cpp")
  demp <- readRDS(test_path("testdata/demographic_projection_object_adult.rds"))
  parameters <- readRDS(test_path("testdata/projection_parameters_adult.rds"))

  # apparently this one is not a list, maybe we should fix that
  parameters$incidinput <- as.list(parameters$incidinput)

  test_tmb(demp, parameters, 1970L, 0L, run_child_model = FALSE)
})
