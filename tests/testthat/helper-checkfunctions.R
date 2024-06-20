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