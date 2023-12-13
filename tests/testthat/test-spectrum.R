test_that("DemProj only matches EPP-ASM", {

  ## Check that population age 15:79 matches between
  ## Note: the open 80+ population does not match because EPP-ASM did
  ##   not handle survivorship of the open age group correctly. This
  ##   is corrected in leapfrog.

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp2 <- prepare_leapfrog_projp(pjnz1)

  lmod1 <- leapfrogR(demp1, hivp1)

  expect_warning(fp1 <- eppasm::prepare_directincid(pjnz1),
                 "no non-missing arguments to min; returning Inf")
  fp1$tARTstart <- 61L

  mod1 <- eppasm::simmod(fp1)

  expect_equal(lmod1$totpop1[16:80,,], mod1[1:65,,1,])

  ## v6.28 -- net-migration at end year
  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_demproj-only_spectrum-v6.28_2023-12-12.PJNZ")
  
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)

  lmod2 <- leapfrogR(demp2, hivp2)

  expect_warning(fp2 <- eppasm::prepare_directincid(pjnz2),
                 "no non-missing arguments to min; returning Inf")
  fp2$tARTstart <- 62L

  mod2 <- eppasm::simmod(fp2)

  expect_equal(lmod2$totpop1[16:80,,], mod2[1:65,,1,])
  
})

test_that("Leapfrog matches DemProj projection without migration", {

  ## v6.13 -- net-migration half at start / half end year
  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff1 <- lmod1$totpop[,,2] - demp1$basepop[,,2]

  specres1 <- eppasm::read_hivproj_output(pjnz1)

  expect_true(all(abs(diff1 < 0.001)))

  ## deaths by sex/age
  expect_true(all(abs(lmod1$natdeaths[,,-1] - specres1$natdeaths[,,-1]) < 0.003))

  ## births by age
  expect_true(all(abs(lmod1$births[-1] - specres1$births[-1]) < 0.002))


  ## v6.28 -- net-migration at end year
  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_demproj-only-no-mig_spectrum-v6.28_2023-12-12.PJNZ")
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)
  lmod2 <- leapfrogR(demp2, hivp2)

  diff2 <- lmod2$totpop[,,2] - demp2$basepop[,,2]

  specres2 <- eppasm::read_hivproj_output(pjnz2)

  expect_true(all(abs(diff2 < 0.001)))

  ## deaths by sex/age
  expect_true(all(abs(lmod2$natdeaths[,,-1] - specres2$natdeaths[,,-1]) < 0.003))

  ## births by age
  expect_true(all(abs(lmod2$births[-1] - specres2$births[-1]) < 0.002))
  
})

test_that("Leapfrog matches DemProj projection with migration", {

  ## v6.13 -- net-migration half at start / half end year
  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff1 <- lmod1$totpop1[,,2:6] - demp1$basepop[,,2:6]

  expect_true(all(abs(diff1) < 0.01))

  ## v6.28 -- net-migration half at start / half end year
  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_demproj-only_spectrum-v6.28_2023-12-12.PJNZ")
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)
  lmod2 <- leapfrogR(demp2, hivp2)

  diff2 <- lmod2$totpop1[,,2:6] - demp2$basepop[,,2:6]

  expect_true(all(abs(diff2) < 0.01))
})
