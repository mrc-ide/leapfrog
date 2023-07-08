test_that("DemProj only matches EPP-ASM", {

  ## Check that population age 15:79 matches between
  ## Note: the open 80+ population does not match because EPP-ASM did
  ##   not handle survivorship of the open age group correctly. This
  ##   is corrected in leapfrog.

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_demproj-only_spectrum-v6.18_2023-07-08.PJNZ")
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  
  lmod <- leapfrogR(demp, hivp)

  expect_warning(fp <- eppasm::prepare_directincid(pjnz1),
                 "no non-missing arguments to min; returning Inf")
  fp$tARTstart <- 61L
    
  mod <- eppasm::simmod(fp)

  expect_equal(lmod$totpop1[16:80,,], mod[1:65,,1,])
})

test_that("Leapfrog matches DemProj projection without migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_demproj-only-no-mig_spectrum-v6.18_2023-07-08.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop[,,2] - demp1$basepop[,,2]

  specres <- eppasm::read_hivproj_output(pjnz1)

  expect_true(all(abs(diff < 0.001)))

  ## deaths by sex/age
  expect_true(all(abs(lmod1$natdeaths[,,-1] - specres$natdeaths[,,-1]) < 0.003))

  ## births by age
  expect_true(all(abs(lmod1$births[-1] - specres$births[-1]) < 0.002))
  
})

test_that("Leapfrog matches DemProj projection with migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop1[,,2:6] - demp1$basepop[,,2:6]

  expect_true(all(abs(diff) < 0.01))
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

  skip("HIV projection with migration is currently failing!")
  
  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_aim-adult-no-art_spectrum-v6.18_2023-07-08.PJNZ")
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
