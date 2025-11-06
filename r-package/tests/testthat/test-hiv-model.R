test_that("stratified HIV population equals single-age HIV population", {

  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  out_full <- run_model(parameters, "HivFullAgeStratification")
  out_coarse <- run_model(parameters, "HivCoarseAgeStratification") 
  out_child <- run_model(parameters, "ChildModel")

  expect_equal(colSums(out_full$p_hivpop[16:81,,]),
               colSums(out_full$h_hivpop,,2) + colSums(out_full$h_artpop,,3))

  expect_equal(colSums(out_coarse$p_hivpop[16:81,,]),
               colSums(out_coarse$h_hivpop,,2) + colSums(out_coarse$h_artpop,,3))

  expect_equal(colSums(out_child$p_hivpop[16:81,,]),
               colSums(out_child$h_hivpop,,2) + colSums(out_child$h_artpop,,3))

})
                 
