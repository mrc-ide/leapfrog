test_that("DemProj only matches EPP-ASM", {

  ## Check that population age 15:79 matches between
  ## Note: the open 80+ population does not match because EPP-ASM did
  ##   not handle survivorship of the open age group correctly. This
  ##   is corrected in leapfrog.

  pjnz1 <- test_path("../testdata/spectrum/v6.29/bwa_demproj-only_spectrum-v6.29_2023-07-08.PJNZ")
  
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

  pjnz1 <- test_path("../testdata/spectrum/v6.29/bwa_demproj-only-no-mig_spectrum-v6.29_2023-07-08.PJNZ")
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

  pjnz1 <- test_path("../testdata/spectrum/v6.29/bwa_demproj-only_spectrum-v6.29_2023-07-08.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop1[,,2:6] - demp1$basepop[,,2:6]

  expect_true(all(abs(diff) < 0.01))
})


test_that("Leapfrog HIV simulation matches EPP-ASM, no ART & no migration", {

  pjnz <- test_path("../testdata/spectrum/v6.29/bwa_aim-adult-no-art_no-migration_spectrum-v6.29_2023-07-08.PJNZ")
  
  demp <- prepare_leapfrog_demp(pjnz)
  hivp <- prepare_leapfrog_projp(pjnz)
  
  lmodC <- leapfrogR(demp, hivp, hiv_strat = "coarse")

  expect_warning(fp <- eppasm::prepare_directincid(pjnz),
                 "no non-missing arguments to min; returning Inf")
  fp$tARTstart <- 61L
    
  mod <- eppasm::simmod(fp)

  ## Checking for age 15-49; 50+ will be slightly distorted by the 80+ projection discrepancy
  expect_equal(lmodC$totpop1[16:50,,], apply(mod[1:35,,,], c(1, 2, 4), sum))
  expect_equal(lmodC$hivpop1[16:50,,], mod[1:35,,2,])
  expect_equal(lmodC$infections[16:81,,], attr(mod, "infections"))
})


test_that("Leapfrog matches AIM projection with no ART and no migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.29/bwa_aim-adult-no-art_no-migration_spectrum-v6.29_2023-07-08.PJNZ")
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


test_that("Leapfrog HIV simulation matches EPP-ASM, no ART & WITH migration", {
  
  pjnz <- test_path("../testdata/spectrum/v6.29/bwa_aim-adult-no-art_spectrum-v6.29_2023-07-08.PJNZ")
  
  demp <- prepare_leapfrog_demp(pjnz)
  hivp <- prepare_leapfrog_projp(pjnz)
  
  expect_warning(fp <- eppasm::prepare_directincid(pjnz),
                 "no non-missing arguments to min; returning Inf")
  fp$tARTstart <- 61L

  ## EPP-ASM is doing something odd with prevalence among child net migrants
  ## Work around this by setting child migration to zero for EPP-ASM validation
  demp$netmigr[1:15,,] <- 0.0
  fp$cumnetmigr[] <- 0.0

  lmodC <- leapfrogR(demp, hivp, hiv_strat = "coarse")
  mod <- eppasm::simmod(fp)

  ## Checking for age 15-49; 50+ will be slightly distorted by the 80+ projection discrepancy
  expect_equal(lmodC$totpop1[16:50,,], apply(mod[1:35,,,], c(1, 2, 4), sum))
  expect_equal(lmodC$hivpop1[16:50,,], mod[1:35,,2,])
  expect_equal(lmodC$infections[16:81,,], attr(mod, "infections"), tolerance = 1e-7)
})        

test_that("Leapfrog matches AIM projection with no ART and WITH migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_aim-adult-no-art_spectrum-v6.18_2023-07-08.PJNZ")
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
  expect_true(all(abs(diff_natdeaths[1:80,,2:61]) < 0.0015))
  expect_true(all(abs(diff_natdeaths[81,,2:61]) < 0.1))

  diff_hivdeaths <- lmod1$hivdeaths - specres$hivdeaths
  expect_true(all(abs(diff_hivdeaths[1:80,,2:61]) < 0.001))
  expect_true(all(abs(diff_hivdeaths[81,,2:61]) < 0.15))
  
})
