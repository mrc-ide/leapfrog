test_that("Leapfrog matches single year age group and coarse age group projection without migration", {

  pjnz1 <- "../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 0.05, threshold_births = 1e-3)
  demog_matches_totpop(pjnz1)
  matches_coarse_age_groups(pjnz1, threshold_pid = c(0, 0, 0), threshold_naturaldeaths = 1e-3)

})

test_that("v6.28 works..", {
    pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")

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
  fp2$tARTstart <- 62L

  mod2 <- eppasm::simmod(fp2)

  expect_equal(lmod2$totpop1[16:80,,], mod2[1:65,,1,])

})

test_that("v6.28 works..", {
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

test_that("v6.28 works..", {

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

##Currently have the right number of births to HIV+ women (checking that with output from spectrum variable hivpregwomen)
test_that('Perinatal transmission of HIV', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_MTCT_perinatal.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp <- prepare_hc_leapfrog_projp(pjnz1, hivp)
  hivp <- leapfrog_input_mods(hivp)
  hivp$ctx_effect <- 0

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  ##look to make sure understanding about HIV incident infections
  ##fertility discounting differences, change all of the FRR to one and see if it lines up
  ##should then just reflect prevalence by age
  ##then account for age, then CD4
  ##check that asfr is divided by 5 in the spectrum inputs file


  lmod <- leapfrogR(demp, hivp)
  specres <- read_hivproj_output(pjnz1)

  x <- data.table(lmod = as.vector(lmod$hiv_births), spec = specres$hivpregwomen)
  x[,diff := spec - lmod]
  x[,ratio := spec/  lmod]
  x[,year := 1970:2030]

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_MTCT_perinatal_pop1.xlsx", ages =0:14, country = 'Botswana')
  #y <- df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_MTCT_perinatal_pop1.xlsx", ages =15, country = 'Botswana')

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

test_that('Perinatal transmission of HIV, some pmtct', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_MTCT_perinatal_options.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0

  hivp$artpaeds_isperc[] <- TRUE

  ##I have no idea what these are
  ## hivp$scalar_art[] <- 1
  ##hivp$fert_rat[] <- 1

  hivp$pmtct[] <- 0
  hivp$pmtct[1:4,which(1970:2030 %in% c(1980:1991)),2] <- 0.1
  hivp$pmtct_mtct[,,2] <- 0

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  hivp$paed_cd4_dist <- c(0.515952304221721, 0.159523042217209, 0.114405414115372, 0.088623912342894, 0.0615533354817918, 0.0380277151144054, 0.0219142765066065)

  ##look to make sure understanding about HIV incident infections
  ##fertility discounting differences, change all of the FRR to one and see if it lines up
  ##should then just reflect prevalence by age
  ##then account for age, then CD4
  ##check that asfr is divided by 5 in the spectrum inputs file

  ## turn off hiv mort
  hivp$paed_cd4_mort[] <- 0
  ## turn off ART
  hivp$paed_art_val[] <- 0
  hivp$adol_cd4_mort[] <- 0
  hivp$cd4_mort_full[] <- 0
  hivp$cd4_mort_coarse[] <- 0

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_MTCT_perinatal_options_pop1.xlsx", ages =0:14, country = 'Botswana')
  ## df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_MTCT_perinatal_pop1.xlsx", ages =15:49, country = 'Botswana')



  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  ##filtering to earlier than 1991 to remove effects of people having kids who were born with hiv
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 1991)
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  # dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  # dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  # diff_art <- abs(dplyr::select(dt_onart, diff))
  # expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
  #


})


test_that('BF transmission of HIV', {
  pjnz <- "../testdata/spectrum/v6.13/TEST_MTCT_bf.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0

  hivp$artpaeds_isperc[] <- TRUE

  ##I have no idea what these are
  ## hivp$scalar_art[] <- 1
  ##hivp$fert_rat[] <- 1

  hivp$pmtct[] <- 0
  hivp$pmtct_mtct[,,1] <- 0
  hivp$pmtct_mtct[,2:5,2] <- 0

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  hivp$paed_cd4_dist <- c(0.515952304221721, 0.159523042217209, 0.114405414115372, 0.088623912342894, 0.0615533354817918, 0.0380277151144054, 0.0219142765066065)

  ##look to make sure understanding about HIV incident infections
  ##fertility discounting differences, change all of the FRR to one and see if it lines up
  ##should then just reflect prevalence by age
  ##then account for age, then CD4
  ##check that asfr is divided by 5 in the spectrum inputs file

  ## turn off hiv mort
  hivp$paed_cd4_mort[] <- 0
  ## turn off ART
  hivp$paed_art_val[] <- 0
  hivp$adol_cd4_mort[] <- 0
  hivp$cd4_mort_full[] <- 0
  hivp$cd4_mort_coarse[] <- 0

  lmod <- leapfrogR(demp, hivp)

  lmod_out <- lmod_output_paed(lmod = lmod)
  df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_MTCT_bf_pop1.xlsx", ages =0:14, country = 'Botswana')
  ## df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_MTCT_perinatal_pop1.xlsx", ages =15:49, country = 'Botswana')



  dt <- dplyr::left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% dplyr::filter(!is.na(pop)) %>% unique()
  ##filtering to earlier than 1991 to remove effects of people having kids who were born with hiv
  dt <- dt %>% dplyr::mutate(diff = lfrog - pop) %>% unique() %>% filter(year < 1991)
  diff = dt$diff
  expect_true(all(abs(diff) < 1e-2), label = 'Off treatment paediatric population in leapfrog and spectrum match')


  # dt_onart <- dplyr::left_join(lmod_out$art, df_out$on_treatment)
  # dt_onart <- dt_onart%>% dplyr::filter(!is.na(pop)) %>% dplyr::mutate(diff = lfrog - pop) %>%  dplyr::ungroup()
  # diff_art <- abs(dplyr::select(dt_onart, diff))
  # expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
  #


})


