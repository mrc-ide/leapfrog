test_that("model runs", {
  pjnz1 <- test_path("testdata/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")

  demp <- leapfrog::prepare_leapfrog_demp(pjnz1)
  out <- run_base_model(demp, 0L)

  ## TODO: expect initial state set up correctly
})