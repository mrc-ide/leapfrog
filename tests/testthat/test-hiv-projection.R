test_that("Adult HIV projection components are consistent", {
  
  pjnz <- test_path("../testdata/spectrum/v6.18/bwa_aim-adult-art-no-special-elig_v6.18_2023-07-08.PJNZ")
  demp <- prepare_leapfrog_demp(pjnz)
  hivp <- prepare_leapfrog_projp(pjnz)
  lmodF <- leapfrogR(demp, hivp, hiv_strat = "full")

  lmodC <- leapfrogR(demp, hivp, hiv_strat = "coarse")

  expect_equal(lmodF$hivpop1[16:81,,],
               colSums(lmodF$hivstrat_adult,,1) + colSums(lmodF$artstrat_adult,,2))

  expect_equal(colSums(lmodC$hivpop1[16:81,,]),
               colSums(lmodC$hivstrat_adult,,2) + colSums(lmodC$artstrat_adult,,3))

  expect_equal(lmodF$hivdeaths[16:81,,],
               colSums(lmodF$aidsdeaths_noart,,1) + colSums(lmodF$aidsdeaths_art,,2))

  expect_equal(colSums(lmodF$hivdeaths[16:81,,]),
               colSums(lmodF$aidsdeaths_noart,,2) + colSums(lmodF$aidsdeaths_art,,3))
  
})
