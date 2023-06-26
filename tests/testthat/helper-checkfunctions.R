#' Checks that coarse age group calculations align between leapfrog and EPP-ASM
#'
#' 
#' @param pjnz The PJNZ file exported from the Spectrum software v6.13
#' @param threshold_pid Absolute threshold of permissible differences between prevalence, incidence, and deaths due to HIV between leapfrog and EPP-ASM. Goal is to minimize these.
#' @param threshold_naturaldeaths Absolute threshold of permissible differences in leapfrog and EPP-ASM's natural deaths. 
#'
#' @return Information on whether models align appropriately
#'

matches_coarse_age_groups <- function(pjnz = "../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ", threshold_pid = c(900, 1, 0.05), 
                                      threshold_naturaldeaths = 2){
  pjnz1 <- test_path(pjnz)
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp, hiv_strat = "coarse")
  
  expect_warning(fp <- eppasm::prepare_directincid(pjnz1),
                 "no non-missing arguments to min; returning Inf")
  fp$tARTstart <- 61L
  
  ## Replace ASFR because demp$asfr is normalised, but fp$asfr is not
  fp$asfr <- demp$asfr
  
  mod <- eppasm::simmod(fp)
  
  #Not doing the open age group because they are calculated differently 
  expect_true(all(abs(attr(mod, 'hivpop')[,1:8,,-61] - lmod$hivstrat_adult[,1:8,,-61]) <= threshold_pid[1]), label = 'Coarse calculation of HIV population matches EPPASM')
  expect_true(all(abs(attr(mod, 'infections')[-66,,] - lmod$infections[16:80,,]) <= threshold_pid[2]), label = 'Coarse calculation of HIV infections matches EPPASM')
  expect_true(all(abs(attr(mod, 'hivdeaths')[-66,,] - lmod$hivdeaths[16:80,,]) <= threshold_pid[3]), label = 'Coarse calculation of HIV deaths matches EPPASM')
  expect_true(all(abs(attr(mod, 'natdeaths')[-66,,] - lmod$natdeaths[16:80,,]) <= threshold_naturaldeaths), label = 'Coarse calculation of non-HIV deaths matches EPPASM')


  
}

demog_matches_totpop <- function(pjnz, threshold = 0.01){
  pjnz1 <- test_path(pjnz)
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)
  
  diff <- lmod1$totpop1 - demp1$basepop
  
  expect_true(all(abs(diff) == 0), label = "Total population and base population align")


  
}

demog_matches_birthsdeaths <- function(pjnz, threshold_deaths = 3, threshold_births = 1e-3){
  pjnz1 <- test_path(pjnz)
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  specres <- eppasm::read_hivproj_output(pjnz1)

  
  ## deaths by sex/age
  expect_true(all(abs(lmod1$natdeaths[,,-1] - specres$natdeaths[,,-1]) < threshold_deaths), 
              label = paste0("Spectrum and leapfrog natural deaths differ by less than ", threshold_deaths, " for all age, sex, year combinations"))

  ## births by age, changed to pct diff
  expect_true(all(abs(lmod1$births[-1] - specres$births[-1]) < threshold_births),
              label = paste0("Spectrum and leapfrog births differ by less than ", threshold_births , " for all years"))
  

}

transmission_matches <- function(pjnz, threshold_absolute_pid = c(250, 25, 3)){
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- test_path(pjnz)
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  
  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)

  specres <- eppasm::read_hivproj_output(pjnz1)

  ## PREVALENCE
  expect_true(all(abs(lmod$hivpop1[1:81,,-1] - specres$hivpop[1:81,,-1]) < threshold_absolute_pid[1]),
              label = paste0("HIV population differs by less than ", threshold_absolute_pid[1], " for all 15+, both sexes, and all years."))

  ##INCIDENCE
  expect_true(all(abs(lmod$infections[1:81,,-1] - specres$infections[1:81,,-1]) < threshold_absolute_pid[2]),
              label = paste0("New infections differ by less than ", threshold_absolute_pid[2], " between leapfrog and Spectrum for 15+, both sexes, and all years."))
  
  ##HIV DEATHS
  expect_true(all(abs(lmod$hivdeaths[1:81,,-1] - specres$hivdeaths[1:81,,-1]) < threshold_absolute_pid[3]),
              label = paste0("HIV deaths in leapfrog differ by less than ", threshold_absolute_pid[3], " from Spectrum for 15+, both sexes, and all years."))

}

spectrum_output <- function(file = "../testdata/spectrum/v6.13/bwa_aim-adult-child-input-art-elig_spectrum-v6.13_2022-02-12_pop1.xlsx", ages = 0:14, country = 'Botswana'){
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself
  df <- file
  if(grepl(pattern = 'testdata', file)){
    df <- test_path(df)
  }
  df <- eppasm::read_pop1(df, country, years = 1970:2022)
  if(any(0:14 %in% ages)){
    df_paed <- df %>% dplyr::filter(age < 5) %>%
      dplyr::right_join(y = data.frame(cd4 = 1:8, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'))) %>%
      dplyr::right_join(y = data.frame(artdur = 2:8, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg') 
     
    df_adol <- df %>% dplyr::filter(age > 4 & age < 15) %>%
      dplyr::right_join(y = data.frame(cd4 = 3:8, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200'))) %>%
      dplyr::right_join(y = data.frame(artdur = 2:8, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg') 
    
    df <- rbind(df_paed, df_adol)
    df <- df %>% dplyr::filter(age %in% ages)
  }else{
    df <- df %>% dplyr::filter(age %in%  c(15:max(ages))) %>%
      dplyr::right_join(y = data.frame(cd4 = 2:8, cd4_cat = c('gte500', '350-500', '250-349', '200-249', '100-199','50-99', 'lte50'))) %>%
      dplyr::right_join(y = data.frame(artdur = 2:8, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg') 
  }

  df <- df %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission)
  
  df_on_treatment <- df[grepl('ART', df$transmission),] %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission) %>% dplyr::group_by(sex, age, cd4_cat, transmission, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique() %>% dplyr::filter(age %in% ages)
  df_off_treatment <- df[!grepl('ART', df$transmission),]%>% dplyr::filter(age %in% ages)
  df_total <- df %>% dplyr::select(sex, age, cd4_cat, year, pop) %>% dplyr::group_by(sex, age, cd4_cat, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique()%>% dplyr::filter(age %in% ages)
  
  return(list(on_treatment = df_on_treatment, off_treatment = df_off_treatment, total = df_total))
  
}

lmod_output_paed <- function(lmod){
  
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
  
  strat_pop_total <- strat_pop %>% dplyr::select(sex, age, cd4_cat, year, lfrog) %>% dplyr::group_by(sex, age, cd4_cat, year) %>% dplyr::mutate(lfrog = sum(lfrog))%>% dplyr::filter(age < 15)
  
  strat_art <- lmod$artstrat_paeds
  dimnames(strat_art) <- list(transmission = c('ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'), cd4_cat =  c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), 
                              age = 0:14, sex = c('Male', 'Female'), year = 1970:2030)
  strat_art <- strat_art %>% as.data.frame.table(responseName = "lfrog")
  strat_art$transmission <- as.character(strat_art$transmission) ; strat_art$cd4_cat <- as.character(strat_art$cd4_cat) ; strat_art$age = as.integer(as.character(strat_art$age)) ; strat_art$sex <- as.character(strat_art$sex) ; strat_art$year <- as.numeric(as.character(strat_art$year))
  strat_art_paed <- strat_art %>% dplyr::filter(age < 5)
  
  strat_art_adol <- strat_art %>% dplyr::filter(age > 4)
  strat_art_adol <- strat_art_adol %>% dplyr::right_join(y = data.frame(cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5'), cd4_cat_new =  c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200', 'lte200')))
  strat_art_adol <-  strat_art_adol %>% dplyr::select(sex, age, cd4_cat = cd4_cat_new, year, lfrog, transmission)
  strat_art_adol <- strat_art_adol %>% dplyr::group_by(sex, age, cd4_cat, year, transmission) %>% dplyr::mutate(lfrog = sum(lfrog)) %>% unique
  strat_art <- rbind(strat_art_paed, strat_art_adol)
  strat_art <- strat_art %>% dplyr::select(sex, age, cd4_cat, year, transmission, lfrog) %>% dplyr::group_by(sex, age, cd4_cat, transmission, year) %>% dplyr::mutate(lfrog = sum(lfrog)) %>% unique() %>% dplyr::filter(age < 15)
  
  
  return(list(prev_strat = strat_pop, prev = strat_pop_total, art = strat_art))
}
