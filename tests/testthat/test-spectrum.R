test_that("DemProj only matches EPP-ASM", {

  ## Check that population age 15:79 matches between
  ## Note: the open 80+ population does not match because EPP-ASM did
  ##   not handle survivorship of the open age group correctly. This
  ##   is corrected in leapfrog.

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_demproj-only_spectrum-v6.18_2023-07-08.PJNZ")
  
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)

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
  fp2$tARTstart <- 61L

  mod2 <- eppasm::simmod(fp2)

  expect_equal(lmod2$totpop1[16:80,,], mod2[1:65,,1,])
  
})

test_that("Leapfrog matches DemProj projection without migration", {

  ## v6.18 -- net-migration half at start / half end year
  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_demproj-only-no-mig_spectrum-v6.18_2023-07-08.PJNZ")
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

  ## v6.18 -- net-migration half at start / half end year
  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_demproj-only_spectrum-v6.18_2023-07-08.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff1 <- lmod1$totpop1 - demp1$basepop

  expect_true(all(abs(diff1) < 0.015))

  ## v6.28 -- net-migration half at start / half end year
  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_demproj-only_spectrum-v6.28_2023-12-12.PJNZ")
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)
  lmod2 <- leapfrogR(demp2, hivp2)

  diff2 <- lmod2$totpop1 - demp2$basepop

  expect_true(all(abs(diff2) < 0.016))
})

test_that("Leapfrog matches AIM projection with no ART and no migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_aim-adult-no-art_no-migration_spectrum-v6.18_2023-07-08.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  specres <- eppasm::read_hivproj_output(pjnz1)

  diff_totpop <- lmod1$totpop1 - specres$totpop
  expect_true(all(abs(diff_totpop[1:80, , ]) < 0.02))
  expect_true(all(abs(diff_totpop[81, , ]) < 0.4))  ## Maybe something to check in 80+ population

  diff_hivpop <- lmod1$hivpop1 - specres$hivpop
  expect_true(all(abs(diff_hivpop[1:80,,]) < 0.01))
  expect_true(all(abs(diff_hivpop[81,,]) < 0.4))

  diff_infections <- lmod1$infections - specres$infections
  expect_true(all(abs(diff_infections) < 0.001))

  diff_natdeaths <- lmod1$natdeaths - specres$natdeaths
  expect_true(all(abs(diff_natdeaths[1:80,,2:61]) < 0.001))
  expect_true(all(abs(diff_natdeaths[81,,2:61]) < 0.1))

  diff_hivdeaths <- lmod1$hivdeaths - specres$hivdeaths
  expect_true(all(abs(diff_hivdeaths[1:80,,2:61]) < 0.001))
  expect_true(all(abs(diff_hivdeaths[81,,2:61]) < 0.1))
  
})


test_that("Leapfrog matches AIM projection with no ART and WITH migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.28/bwa_aim-adult-no-art_spectrum-v6.28_2023-12-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  specres <- eppasm::read_hivproj_output(pjnz1)

  diff_totpop <- lmod1$totpop1 - specres$totpop
  expect_true(all(abs(diff_totpop[1:80, , ]) < 0.025))
  expect_true(all(abs(diff_totpop[81, , ]) < 0.5))  ## Maybe something to check in 80+ population

  diff_hivpop <- lmod1$hivpop1 - specres$hivpop
  expect_true(all(abs(diff_hivpop[1:80,,]) < 0.01))
  expect_true(all(abs(diff_hivpop[81,,]) < 0.5))

  diff_infections <- lmod1$infections - specres$infections
  expect_true(all(abs(diff_infections) < 0.001))

  diff_natdeaths <- lmod1$natdeaths - specres$natdeaths
  expect_true(all(abs(diff_natdeaths[1:80,,2:61]) < 0.0011))
  expect_true(all(abs(diff_natdeaths[81,,2:61]) < 0.1))

  diff_hivdeaths <- lmod1$hivdeaths - specres$hivdeaths
  expect_true(all(abs(diff_hivdeaths[1:80,,2:61]) < 0.001))
  expect_true(all(abs(diff_hivdeaths[81,,2:61]) < 0.11))
  
})
