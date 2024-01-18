test_that("Leapfrog matches single year age group and coarse age group projection without migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_demproj-only_spectrum-v6.18_2023-07-08.PJNZ")

  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  hivp1 <- prepare_hc_leapfrog_projp(pjnz1, hivp1)


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
  hivp2 <- prepare_hc_leapfrog_projp(pjnz2, hivp2)

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
  hivp1 <- prepare_hc_leapfrog_projp(pjnz1, hivp1)
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
  hivp2 <- prepare_hc_leapfrog_projp(pjnz2, hivp2)
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
  hivp1 <- prepare_hc_leapfrog_projp(pjnz1, hivp1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff1 <- lmod1$totpop1 - demp1$basepop

  expect_true(all(abs(diff1) < 0.015))

  ## v6.28 -- net-migration half at start / half end year
  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_demproj-only_spectrum-v6.28_2023-12-12.PJNZ")
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)
  hivp2 <- prepare_hc_leapfrog_projp(pjnz2, hivp2)
  lmod2 <- leapfrogR(demp2, hivp2)

  diff2 <- lmod2$totpop1 - demp2$basepop

  expect_true(all(abs(diff2) < 0.016))
})


test_that("Leapfrog matches AIM projection with no ART and no migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.18/bwa_aim-adult-no-art_no-migration_spectrum-v6.18_2023-07-08.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  hivp1 <- prepare_hc_leapfrog_projp(pjnz1, hivp1)
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

  ##Maggie added this 08/01/2024
  pjnz2 <- test_path("../testdata/spectrum/v6.28/bwa_aim-adult-no-art_no-migration_spectrum-v6.28_2024-08-01.PJNZ")
  demp2 <- prepare_leapfrog_demp(pjnz2)
  hivp2 <- prepare_leapfrog_projp(pjnz2)
  hivp2 <- prepare_hc_leapfrog_projp(pjnz2, hivp2)
  lmod2 <- leapfrogR(demp2, hivp2)

  specres2 <- eppasm::read_hivproj_output(pjnz2)


})


test_that("Leapfrog matches AIM projection with no ART and WITH migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.28/bwa_aim-adult-no-art_spectrum-v6.28_2023-12-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  hivp1 <- prepare_hc_leapfrog_projp(pjnz1, hivp1)
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

##TODO: add in test for hiv entrant population
# test_that("Age 15 entrant population matches", {
#
#
# })

test_that("Leapfrog matches direct incidence option at coarse and single year age structure, no ART and no hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-no-hiv-deaths_spectrum-v6.13_2022-02-12.pjnz"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-1, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(0.2, 1e-3, 1e-3))
  #matches_coarse_age_groups(pjnz1, threshold_pid = c(0.5, 0.01, 1e-3))
})

test_that("Leapfrog matches direct incidence option at coarse and single year age structure, no ART + hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-1, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(0.5, 1e-3, 0.1))
#  matches_coarse_age_groups(pjnz1, threshold_pid = c(27, 4, 3), threshold_naturaldeaths = 3)

})

test_that("Input childhood infections and test alignment betwen leapfrog and spectrum", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models

  ##100 children under 5 getting infected in 1980
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 3e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(0.5, 1e-3, 1e-3))
##  matches_coarse_age_groups(pjnz1, threshold_pid = c(11.5, 0.2, 1e-3))
})

test_that('Paediatric aging and natural deaths working appropriately', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- spectrum_output(file = df, ages = 0:14, country = 'botswana')
  df <- data.frame(df$off_treatment)

  ##5-15 can move more than one cd4 category in a year... feels wrong?
  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  strat_pop_paed <- strat_pop %>% dplyr::filter(age < 5)

  strat_pop_adol <- strat_pop %>% dplyr::filter(age > 4)
  strat_pop_adol <- strat_pop_adol %>% dplyr::right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_pop_adol <-  strat_pop_adol %>% dplyr::select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_pop_adol <- strat_pop_adol %>% dplyr::group_by(sex, age, cd4_cat, year, transmission) %>% dplyr::mutate(lfrog = sum(lfrog)) %>% unique
  strat_pop <- rbind(strat_pop_paed, strat_pop_adol)

  dt <- dplyr::right_join(df, strat_pop)
  dt <- dt %>% dplyr::filter(!is.na(pop))
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop)

  ##check that the populations between specrum and lfrog match
  expect_true(all(dplyr::select(dt, diff) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  })

test_that('Paediatric transition through CD4 working appropriately', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-transition_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-transition_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- spectrum_output(file = df, ages = 0:14, country = 'botswana')
  df <- data.frame(df$off_treatment)

  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-14', '5-10', 'lte5'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  strat_pop_paed <- strat_pop %>% dplyr::filter(age < 5)

  strat_pop_adol <- strat_pop %>% dplyr::filter(age > 4)
  strat_pop_adol <- strat_pop_adol %>% dplyr::right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-14', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_pop_adol <-  strat_pop_adol %>% dplyr::select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_pop_adol <- strat_pop_adol %>% dplyr::group_by(sex, age, cd4_cat, year, transmission) %>% dplyr::mutate(lfrog = sum(lfrog)) %>% unique
  strat_pop <- rbind(strat_pop_paed, strat_pop_adol)

  dt <- dplyr::right_join(df, strat_pop)
  dt <- dt %>% dplyr::filter(!is.na(pop))

  dt <- dt %>% dplyr::mutate(diff = lfrog - pop)
  expect_true(all(abs(dplyr::select(dt, diff)) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')


})

test_that('Paediatric HIV mortality working as expected', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-hivmort_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)


  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-hivmort_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- spectrum_output(file = df, ages = 0:14, country = 'botswana')
  df <- data.frame(df$off_treatment)

  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  strat_pop_paed <- strat_pop %>% dplyr::filter(age < 5)

  strat_pop_adol <- strat_pop %>% dplyr::filter(age > 4)
  strat_pop_adol <- strat_pop_adol %>% dplyr::right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_pop_adol <-  strat_pop_adol %>% dplyr::select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_pop_adol <- strat_pop_adol %>% dplyr::group_by(sex, age, cd4_cat, year, transmission) %>% dplyr::mutate(lfrog = sum(lfrog)) %>% unique
  strat_pop <- rbind(strat_pop_paed, strat_pop_adol)

  dt <- dplyr::right_join(df, strat_pop)
  dt <- dt %>% dplyr::filter(!is.na(pop))

  ## 13 and 14 year old females are causing NANs are
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop)
  expect_true(all(abs(dplyr::select(dt, diff)) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')


})

test_that('ART % implemented, mortality reduction & all eligible', {
  pjnz <- "../testdata/spectrum/v6.13/TEST.PJNZ"
  pjnz1 <- test_path(pjnz)
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))
  hivp$paed_art_val <- hivp$paed_art_val / 100

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')



})

test_that('ART counts implemented, mortality reduction & all eligible', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_art_COUNTS.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))
  hivp$artpaeds_isperc[] <- FALSE


  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_art_COUNTS_pop1.xlsx", ages =0:14, country = 'Botswana')



  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')


})

test_that('ART counts, number covered is less than total prevalent cases', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_art_COUNTS_insufficient.PJNZ"
  pjnz1 <- test_path(pjnz)
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))
  hivp$artpaeds_isperc[] <- FALSE

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_art_COUNTS_insufficient_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')


})

test_that('ART counts, number covered is less than total prevalent cases. ART coverage goes down', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_art_COUNTS_decreasing_ART.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))
  hivp$artpaeds_isperc[] <- FALSE

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_art_COUNTS_decreasing_ART_pop1.xlsx", ages =0:14, country = 'Botswana')


  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
})

test_that('ART counts, number covered is less than total prevalent cases. ART input switches from numbers to percentages', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_art_COUNTS_num_to_pct.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  hivp$paed_art_val[which(1970:2030 %in% 1997:2002)] <- hivp$paed_art_val[which(1970:2030 %in% 1997:2002)]/100
  hivp$artpaeds_isperc[] <- FALSE
  hivp$artpaeds_isperc[which(1970:2030 %in% 1997:2002)] <- TRUE
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_art_COUNTS_num_to_pct_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')

  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(abs(diff_art) < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
})

test_that('BWA normal treatment', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_art_BWA.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)

  hivp$ctx_effect <- 0
  ##quick fix for nosocomial infections
  hivp$paed_cd4_dist <- c(1,rep(0,6))
  hivp$paed_art_val <- hivp$paed_art_val / 100


  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_art_BWA_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')

  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
})

##08/01/2024 working
test_that('Perinatal transmission of HIV', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_perinatal.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$ctx_effect <- 0
  hivp$paed_cd4_transition[6,1] <- 1 - sum(hivp$paed_cd4_transition[-6,1])
  hivp$paed_cd4_transition[6,7] <- 1 - sum(hivp$paed_cd4_transition[-6,7])
  hivp$adult_cd4_dist[3,5] <- 0.6665589
  hivp$adult_cd4_dist[4,5] <- 1 - hivp$adult_cd4_dist[3,5]


  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)


## issue is that lmod doesn't have any kids? in hivpop1?
  lmod <- leapfrogR(demp, hivp)
  specres <- read_hivproj_output(pjnz1)

  x <- data.table(lmod = as.vector(lmod$hiv_births), spec = specres$hivpregwomen)
  x[,diff := spec - lmod]
  x[,ratio := spec/  lmod]
  x[,year := 1970:2030]

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_perinatal_pop1.xlsx", ages =0:14, country = 'Botswana')
  y <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_perinatal_pop1.xlsx", ages =15, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique()

  x <- data.table(dt)
  x[year == 2000 & age == 14 & sex == 'Female']
  y <- data.table(y$total)
  y <- y[year == 2001 & age == 15 & sex == 'Female']
  y
  z <- lmod$hivstrat_adult[,1,2,which(1970:2030 == 2001)]
  z

  diff = dt$diff[!is.nan(dt$diff)]
  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  diff_art <- abs(dplyr::select(dt_onart, diff))
  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
})

##09/01/2024 working
test_that('Perinatal transmission of HIV, some pmtct', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_perinatal_options.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$ctx_effect <- 0
  hivp$paed_cd4_transition[6,1] <- 1 - sum(hivp$paed_cd4_transition[-6,1])
  hivp$paed_cd4_transition[6,7] <- 1 - sum(hivp$paed_cd4_transition[-6,7])
  hivp$adult_cd4_dist[3,5] <- 0.6665589
  hivp$adult_cd4_dist[4,5] <- 1 - hivp$adult_cd4_dist[3,5]
  hivp$pmtct_input_isperc[] <- T

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)


  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  specres <- read_hivproj_output(pjnz1)
  x <- data.table(lmod = as.vector(lmod$hiv_births), spec = specres$hivpregwomen)
  x[,diff := spec - lmod]
  x[,ratio := spec/  lmod]
  x[,year := 1970:2030]

  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_perinatal_options_pop1.xlsx", ages =0:14, country = 'Botswana')
  ## y <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_perinatal_pop1.xlsx", ages =15, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')


})

##09/01/2024 working
test_that('BF transmission of HIV', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_bf.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$ctx_effect <- 0
  hivp$paed_cd4_transition[6,1] <- 1 - sum(hivp$paed_cd4_transition[-6,1])
  hivp$paed_cd4_transition[6,7] <- 1 - sum(hivp$paed_cd4_transition[-6,7])
  hivp$adult_cd4_dist[3,5] <- 0.6665589
  hivp$adult_cd4_dist[4,5] <- 1 - hivp$adult_cd4_dist[3,5]
  hivp$pmtct_input_isperc[] <- T


  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  specres <- read_hivproj_output(pjnz1)
  x <- data.table(lmod = as.vector(lmod$hiv_births), spec = specres$hivpregwomen)
  x[,diff := spec - lmod]
  x[,ratio := spec/  lmod]
  x[,year := 1970:2030]

  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_bf_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev_strat, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')

})

##09/01/2024 working
test_that('BF transmission of HIV with pmtct', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_bf_options.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$ctx_effect <- 0
  hivp$paed_cd4_transition[6,1] <- 1 - sum(hivp$paed_cd4_transition[-6,1])
  hivp$paed_cd4_transition[6,7] <- 1 - sum(hivp$paed_cd4_transition[-6,7])
  hivp$adult_cd4_dist[3,5] <- 0.6665589
  hivp$adult_cd4_dist[4,5] <- 1 - hivp$adult_cd4_dist[3,5]
  hivp$pmtct_input_isperc[] <- T


  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  specres <- read_hivproj_output(pjnz1)
  x <- data.table(lmod = as.vector(lmod$hiv_births), spec = specres$hivpregwomen)
  x[,diff := spec - lmod]
  x[,ratio := spec/  lmod]
  x[,year := 1970:2030]

  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_bf_options_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev_strat, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')

})

##09/01/2024 working
test_that('BF and perinatal transmission of HIV, no pmtct', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$ctx_effect <- 0
  hivp$paed_cd4_transition[6,1] <- 1 - sum(hivp$paed_cd4_transition[-6,1])
  hivp$paed_cd4_transition[6,7] <- 1 - sum(hivp$paed_cd4_transition[-6,7])
  hivp$adult_cd4_dist[3,5] <- 0.6665589
  hivp$adult_cd4_dist[4,5] <- 1 - hivp$adult_cd4_dist[3,5]
  hivp$pmtct_input_isperc[] <- T


  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  specres <- read_hivproj_output(pjnz1)
  x <- data.table(lmod = as.vector(lmod$hiv_births), spec = specres$hivpregwomen)
  x[,diff := spec - lmod]
  x[,ratio := spec/  lmod]
  x[,year := 1970:2030]

  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pop1.xlsx", ages =0:14, country = 'Botswana')

  dt <- dplyr::left_join(lmod_out$prev_strat, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')

})

##17/01/2024 Working
test_that('BF and perinatal transmission of HIV, pmtct', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pmtct.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$art_mtct[,,2] <- hivp$art_mtct[,c(3,2,1),2]

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pmtct_pop1.xlsx", ages =0:14, country = 'Botswana', years_in = 1990:2020)

  dt <- dplyr::left_join(lmod_out$prev_strat, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)

  diff = dt$diff
  expect_true(all(abs(diff) < 1.5e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')

})

##17/01/2024 Working, with a highish threshold for differences
test_that('BF and perinatal transmission of HIV, pmtct, paed mort', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pmtct_paed_mort.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$art_mtct[,,2] <- hivp$art_mtct[,c(3,2,1),2]

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  specres <- read_hivproj_output(pjnz1)


  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pmtct_paed_mort_pop1.xlsx", ages =0:14, country = 'Botswana', years_in = 1970:2020)

  dt <- dplyr::left_join(lmod_out$prev_strat, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)
  x = data.table(dt)
  diff = dt$diff
  ##Note that the higher threhold for differences here is due to the misalignment in total
  ##population causing incidence to be different
  expect_true(all(abs(diff) < 6.5e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')

})

##17/01/2024 Working with no cd4 transition
test_that('BF and perinatal transmission of HIV, pmtct & ctx', {
  pjnz <- "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pmtct_paed_mort_cotrim.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp$art_mtct[,,2] <- hivp$art_mtct[,c(3,2,1),2]

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  lmod_out <- lmod_output_paed(lmod = lmod)

  specres <- read_hivproj_output(pjnz1)

  df_out <- spectrum_output(file = "../testdata/spectrum/v6.28/TEST_MTCT_BF_PERI_pmtct_paed_mort_cotrim_pop1.xlsx", ages =0:14, country = 'Botswana', years_in = 2000:2020)

  dt <- dplyr::left_join(lmod_out$prev_strat, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 2030)
  x = data.table(dt)
  diff = dt$diff
  ##Note that the higher threhold for differences here is due to the misalignment in total
  ##population causing incidence to be different
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')
})







