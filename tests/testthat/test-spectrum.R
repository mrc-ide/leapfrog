test_that("Leapfrog matches single year age group and coarse age group projection without migration", {
  
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 0.05, threshold_births = 1e-3)
  demog_matches_totpop(pjnz1)
  matches_coarse_age_groups(pjnz1, threshold_pid = c(0, 0, 0), threshold_naturaldeaths = 1e-3)
  
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
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(0.2, 1e-3, 1e-3))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(0.5, 0.01, 1e-3))
})

test_that("Leapfrog matches direct incidence option at coarse and single year age structure, no ART + hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(0.5, 1e-3, 0.1))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(27, 4, 3), threshold_naturaldeaths = 3)

})

test_that("Input childhood infections and test alignment betwen leapfrog and spectrum", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  
  ##100 children under 5 getting infected in 1980
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
##  transmission_matches(pjnz1, threshold_absolute_pid = c(0.2, 1e-3, 1e-3))
##  matches_coarse_age_groups(pjnz1, threshold_pid = c(11.5, 0.2, 1e-3))
})


test_that('Paediatric aging and natural deaths working appropriately', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0
  hivp$paed_art_val[] <- 0
  hivp$paed_cd4_mort[] <- 0
  hivp$adol_cd4_mort[] <- 0
  hivp$adol_cd4_prog[] <- 0
  hivp$paed_cd4_prog[] <- 0
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)
  
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- test_path(df)
  df <- read_pop1(df, "Botswana", years = 1970:2022)
  df_paed <- df %>% filter(age < 5) %>%
    right_join(y = data.frame(cd4 = 1:8, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  
  df_adol <- df %>% filter(age > 4) %>%
    right_join(y = data.frame(cd4 = 3:8, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  
  df <- rbind(df_paed, df_adol)
  df <- df %>% select(sex, age, cd4_cat, year, pop, transmission)
  
  
  ##5-15 can move more than one cd4 category in a year... feels wrong?
  
  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  strat_pop_paed <- strat_pop %>% filter(age < 5)
  
  strat_pop_adol <- strat_pop %>% filter(age > 4)
  strat_pop_adol <- strat_pop_adol %>% right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_pop_adol <-  strat_pop_adol %>% select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_pop_adol <- strat_pop_adol %>% group_by(sex, age, cd4_cat, year, transmission) %>% mutate(lfrog = sum(lfrog)) %>% unique
  strat_pop <- rbind(strat_pop_paed, strat_pop_adol)
  
  dt <- right_join(df, strat_pop)
  dt <- dt %>% filter(!is.na(pop))
  
  ## dt <- dt %>% filter(age < 5 & !is.na(pop))
  
  dt <- dt %>% mutate(diff = lfrog - pop)
  
  ##check that the populations between specrum and lfrog match
  expect_true(all(select(dt, diff) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  
  
  })


test_that('Paediatric transition through CD4 working appropriately', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-transition_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0
  hivp$paed_art_val[] <- 0
  hivp$paed_cd4_mort[] <- 0
  hivp$adol_cd4_mort[] <- 0
  
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)
  
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself 
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-transition_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- test_path(df)
  df <- read_pop1(df, "Botswana", years = 1970:2022)
  df_paed <- df %>% filter(age < 5) %>%
    right_join(y = data.frame(cd4 = 1:8, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-14', '5-10', 'lte5'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  
  df_adol <- df %>% filter(age > 4) %>%
    right_join(y = data.frame(cd4 = 3:8, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  
  df <- rbind(df_paed, df_adol)
  df <- df %>% select(sex, age, cd4_cat, year, pop, transmission)
  
  
  ##5-15 can move more than one cd4 category in a year... feels wrong?
  
  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-14', '5-10', 'lte5'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  strat_pop_paed <- strat_pop %>% filter(age < 5)
  
  strat_pop_adol <- strat_pop %>% filter(age > 4)
  strat_pop_adol <- strat_pop_adol %>% right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-14', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_pop_adol <-  strat_pop_adol %>% select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_pop_adol <- strat_pop_adol %>% group_by(sex, age, cd4_cat, year, transmission) %>% mutate(lfrog = sum(lfrog)) %>% unique
  strat_pop <- rbind(strat_pop_paed, strat_pop_adol)
  
  dt <- right_join(df, strat_pop)
  dt <- dt %>% filter(!is.na(pop))
  
  ## dt <- dt %>% filter(age < 5 & !is.na(pop))
  
  dt <- dt %>% mutate(diff = lfrog - pop)
  expect_true(all(abs(select(dt, diff)) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  
})


test_that('Paediatric HIV mortality working as expected', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-hivmort_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0
  hivp$paed_art_val[] <- 0
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)
  
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself 
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-hivmort_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- test_path(df)
  df <- read_pop1(df, "Botswana", years = 1970:2022)
  df_paed <- df %>% filter(age < 5) %>%
    right_join(y = data.frame(cd4 = 1:8, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  
  df_adol <- df %>% filter(age > 4) %>%
    right_join(y = data.frame(cd4 = 3:8, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  
  df <- rbind(df_paed, df_adol)
  df <- df %>% select(sex, age, cd4_cat, year, pop, transmission)
  

  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  strat_pop_paed <- strat_pop %>% filter(age < 5)
  
  strat_pop_adol <- strat_pop %>% filter(age > 4)
  strat_pop_adol <- strat_pop_adol %>% right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_pop_adol <-  strat_pop_adol %>% select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_pop_adol <- strat_pop_adol %>% group_by(sex, age, cd4_cat, year, transmission) %>% mutate(lfrog = sum(lfrog)) %>% unique
  strat_pop <- rbind(strat_pop_paed, strat_pop_adol)
  
  dt <- right_join(df, strat_pop)
  dt <- dt %>% filter(!is.na(pop))
  
  ## 13 and 14 year old females are causing NANs are 
  dt <- dt %>% mutate(diff = lfrog - pop)
  expect_true(all(abs(select(dt, diff)) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  
})



test_that('ART counts implemented, no mortality reduction & all eligible', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-child-input-art-elig_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0
  hivp$paed_art_mort[] <- 0
  hivp$adol_art_mort[] <- 0
  hivp$paed_art_val[which(1970:2030 %in% 1995:2030)] <- 1
  hivp$paed_art_elig_age[] <- 15

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  lmod <- leapfrogR(demp, hivp)
  
  lmod_out <- lmod_output_paed(lmod = lmod)
  x=data.table(lmod_out$prev_strat)
##  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  #df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/bwa_aim-adult-child-input-art-elig_spectrum-v6.13_2022-02-12_pop1.xlsx", ages =0:14, country = 'Botswana')

  

  dt <- left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% filter(!is.na(pop)) %>% unique()
  dt <- dt %>% mutate(diff = lfrog - pop) %>% unique()
 ## diff = dt$diff
  ##expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')
  
  
  dt_onart <- left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% filter(!is.na(pop))  %>% mutate(diff = lfrog - pop) %>%  ungroup()
 ## diff_art <- abs(select(dt_onart, diff))
  ##expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
  x.1 = data.table(dt)
  y.1 = data.table(dt_onart)
  
})

test_that('ART counts implemented, mortality reduction & all eligible', {
  pjnz <- "../testdata/spectrum/v6.13/TEST.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0
  hivp$paed_art_val[which(1970:2030 %in% 1995:2030)] <- 1
  hivp$paed_art_elig_age[] <- 15
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)
  
  lmod_out <- lmod_output_paed(lmod = lmod)
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
 ## df_out <- spectrum_output(file = "../testdata/spectrum/v6.13/TEST_pop1.xlsx", ages =0:14, country = 'Botswana')
  
  
  
  dt <- left_join(lmod_out$prev, df_out$off_treatment)
  dt <- dt %>% filter(!is.na(pop)) %>% unique()
  dt <- dt %>% mutate(diff = lfrog - pop) %>% unique()
  diff = dt$diff
##  expect_true(all(abs(diff) < 1e-3), label = 'Off treatment paediatric population in leapfrog and spectrum match')
  
  
  dt_onart <- left_join(lmod_out$art, df_out$on_treatment)
  dt_onart <- dt_onart%>% filter(!is.na(pop)) %>% mutate(diff = lfrog - pop) %>%  ungroup()
  diff_art <- abs(select(dt_onart, diff))
##  expect_true(all(diff_art < 1e-3), label = 'On treatment paediatric population in leapfrog and spectrum match')
  
  x = data.table(dt)
  y = data.table(dt_onart)
  
})

