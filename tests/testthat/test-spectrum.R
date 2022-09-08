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
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  ##right now paed progress is hardcoded, corrected for this test case here
  hivp$paed_cd4_prog[] <- 0
  hivp$paed_cd4_mort[] <- 0
  
  
  
  lmod <- leapfrogR(demp, hivp)
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- test_path(df)
  df <- read_pop1(df, "Botswana", years = 1970:2022)
  df <- df %>% filter(age < 5) %>%
    right_join(y = data.frame(cd4 = 1:7, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10'))) %>%
    filter(cd4_cat != 'neg') 
  
  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'x'),transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+') , age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- apply(strat_pop, MARGIN = c(1,3:5), FUN = sum)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  
  dt <- right_join(df, strat_pop)
  dt <- dt %>% filter(age < 5 & !is.na(pop))
  dt <- dt %>% mutate(diff = lfrog - pop)
  
  stratified_pop <- apply(lmod$hivstrat_paeds, MARGIN = c(3,5), FUN = sum)
  pop <- apply(lmod$hivpop1[1:15,,], MARGIN = c(1,3), FUN = sum)
  
  
  expect_true(all(stratified_pop - pop < 1e-3), label = 'Stratified population and unstrat pop match for paeds')
  
  ##check that the populations between specrum and lfrog match
  expect_true(all(select(dt, diff) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  
  
  })


test_that('Paediatric transition through CD4 working appropriately', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-transition_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  hivp$paed_cd4_mort[] <- 0
  
  
  lmod <- leapfrogR(demp, hivp)
  
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself 
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-transition_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- test_path(df)
  df <- read_pop1(df, "Botswana", years = 1970:2022)
  df <- df %>% filter(age < 5) %>%
    right_join(y = data.frame(cd4 = 1:7, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10'))) %>%
    filter(cd4_cat != 'neg') 
  df <- df %>% select(sex, age, cd4_cat, year, pop)
  
  
  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'x'),transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+') , age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- apply(strat_pop, MARGIN = c(1,3:5), FUN = sum)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  
  
  dt <- right_join(df, strat_pop)
  dt <- dt %>% filter(age < 5 & !is.na(pop))
  
  dt <- dt %>% mutate(diff = lfrog - pop)
  expect_true(all(select(dt, diff) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  
})


test_that('Paediatric HIV mortality working as expected', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-hivmort_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  hivp$ctx_effect <- 0
  hivp$ctx_val[] <- 0
  
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
  df <- df %>% filter(age < 5) %>%
    right_join(y = data.frame(cd4 = 1:7, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10'))) %>%
    right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
    filter(cd4_cat != 'neg') 
  df <- df %>% select(sex, age, cd4_cat, year, pop, transmission)
  
  
  strat_pop <- lmod$hivstrat_paeds
  dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'x'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
  strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))
  
  
  dt <- right_join(df, strat_pop)
  dt <- dt %>% filter(age < 5 & !is.na(pop))
  
  dt <- dt %>% mutate(diff = lfrog - pop)
  expect_true(all(select(dt, diff) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')
  
})



test_that('Cotrim reduction of mortality implemented', {
  pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-cotrim_spectrum-v6.13_2022-02-12.PJNZ"
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)
  
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself 
  source("https://raw.githubusercontent.com/mrc-ide/eppasm/new-master/R/read-spectrum-pop1.R")
  df <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-cotrim_spectrum-v6.13_2022-02-12_pop1.xlsx"
  df <- test_path(df)
  df <- read_pop1(df, "Botswana", years = 1970:2022)
  df <- df %>% filter(age < 5) %>%
  right_join(y = data.frame(cd4 = 1:7, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10'))) %>%
  right_join(y = data.frame(artdur = 2:5, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'))) %>%
  filter(cd4_cat != 'neg') 
  df <- df %>% select(sex, age, cd4_cat, year, pop, transmission)


strat_pop <- lmod$hivstrat_paeds
dimnames(strat_pop) <- list(cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'x'), transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                            age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
strat_pop <- strat_pop %>% as.data.frame.table(responseName = "lfrog")
strat_pop$cd4_cat <- as.character(strat_pop$cd4_cat) ; strat_pop$age = as.integer(as.character(strat_pop$age)) ; strat_pop$sex <- as.character(strat_pop$sex) ; strat_pop$year <- as.numeric(as.character(strat_pop$year))


dt <- right_join(df, strat_pop)
dt <- dt %>% filter(age < 5 & !is.na(pop))

dt <- dt %>% mutate(diff = lfrog - pop)
expect_true(all(select(dt, diff) < 1e-3), label = 'Prevalence in leapfrog and spectrum match')

})






