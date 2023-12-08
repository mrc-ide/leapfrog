rm(list = ls())
devtools::load_all(".")
library(tidyr)
library(tidyverse)
library(dplyr)
library(data.table)
#library(leapfrog)
devtools::load_all('C:/Users/mwalters/leapfrog')

spectrum_output <- function(file = "../testdata/spectrum/v6.13/bwa_aim-adult-child-input-art-elig_spectrum-v6.13_2022-02-12_pop1.xlsx", ages = 0:14, country = 'Botswana'){
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself
  df <- file

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
      dplyr::right_join(y = data.frame(cd4 = 1:8, cd4_cat = c('hivneg','gte500', '350-500', '250-349', '200-249', '100-199','50-99', 'lte50'))) %>%
      dplyr::right_join(y = data.frame(artdur = 1:8, transmission = c('hivneg','perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg')
  }
  # setnames(df, 'cd4', 'cd4_cat')
  df <- df %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission)

  df_on_treatment <- df[grepl('ART', df$transmission),] %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission) %>% dplyr::group_by(sex, age, cd4_cat, transmission, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique() %>% dplyr::filter(age %in% ages)
  df_off_treatment <- df[!grepl('ART', df$transmission),]%>% dplyr::filter(age %in% ages)
  df_hivpop <- df[df$transmission != 'hivneg',] %>% dplyr::select(sex, age, cd4_cat, year, pop) %>% dplyr::group_by(sex, age, cd4_cat, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique()%>% dplyr::filter(age %in% ages)
  df_total <- df %>% dplyr::select(sex, age, cd4_cat, year, pop) %>% dplyr::group_by(sex, age, cd4_cat, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique()%>% dplyr::filter(age %in% ages)

  return(list(on_treatment = df_on_treatment, off_treatment = df_off_treatment, hivpop = df_hivpop, total = df_total))

}

demp <- readRDS(testthat::test_path("testdata/demographic_projection_object_adult.rds"))
parameters <- readRDS(testthat::test_path("testdata/projection_parameters_adult.rds"))

out <- run_model(demp, parameters, NULL, NULL, 0:60,
                 run_child_model = FALSE, pop_adjust = TRUE)

pjnz1 <- testthat::test_path("testdata/bwa_aim-no-special-elig.PJNZ")
 #parameters <- leapfrog:::prepare_hc_leapfrog_projp(pjnz1, parameters)

 lmod1 <- leapfrogR(demp, parameters)
#
#
#
 saveRDS(lmod1, 'C:/Users/mwalters/Desktop/leapfrog.RDS')
specres <- eppasm::read_hivproj_output(pjnz1)

spec = spectrum_output(file = 'C:/Users/mwalters/frogger/tests/testthat/testdata/bwa_aim-no-special-elig_pop1.xlsx', ages = 15:49, country = 'bwa')

spec_hivpop <- data.table(spec$hivpop)
spec_hivpop <- spec_hivpop[,.(pop = sum(pop)), by = c('year', 'age')]
spec_hivpop <- spec_hivpop[age == 15,]


frog_hivpop <- apply(out$p_hiv_pop[16,,which(1970:2030 %in% spec_hivpop$year)], MARGIN = 2, sum)
spec_hivpop[,frogger := frog_hivpop]
spec_hivpop[,diff := pop - frog_hivpop]
frog_strat_pop <- apply(out$h_hiv_adult[,1,,which(1970:2030 %in% spec_hivpop$year)], MARGIN = c(3), sum)
spec_hivpop[,strat := frog_strat_pop]
spec_hivpop


spec_total <- data.table(spec$total)
spec_total <- spec_total[,.(pop = sum(pop)), by = c('year', 'age', 'sex')]
spec_total <- spec_total[age == 15 & sex == 'Female',]
frog_total <- out$p_total_pop[16,2,which(1970:2030 %in% spec_total$year)]
spec_total[,frogger := frog_total]
spec_total[,diff := pop - frogger]
spec_total

#0.10996644 / 0.10925
#1.006558
demp$targetpop[16,2,2]
