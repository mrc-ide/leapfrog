#!/usr/bin/env Rscript

## This script creates test data required for the frogger tests.
## We read some input data and prepare a set of demographic projection
## and HIV parameters for both the adult and the child model.
## We also run leapfrog and save out the result for use in reference tests

# nolint start
library(frogger)
library(data.table)
library(dplyr)
source('./scripts/spectrum_inputs_paeds.R')
source('./scripts/read_spectrum.R')

## Create demographic and projection parameters for adults
pjnz_adult <- system.file("pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ", package = "frogger", mustWork = TRUE)

# demp <- prepare_leapfrog_demp(pjnz_adult)
# saveRDS(demp, testthat::test_path("testdata/demographic_projection_object_adult.rds"))
# proj <- prepare_leapfrog_projp(pjnz_adult)
# saveRDS(proj, testthat::test_path("testdata/projection_parameters_adult.rds"))
#
# # Used as reference data (Run from leapfrog/master)
# lmod <- leapfrog::leapfrogR(demp, proj)
# saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit.rds"))
#
# lmod <- leapfrog::leapfrogR(demp, proj, hiv_strat = "coarse")
# saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit_coarse.rds"))
#
# mod <- leapfrog::leapfrogR(demp, proj, hiv_steps_per_year = 0L)
# saveRDS(lmod, testthat::test_path("testdata/fit_demography.rds"))


#Create paeds parameters (Run from leapfrog/uncertainrt_analysis_working)
pjnz_child <- testthat::test_path("testdata/bwa_aim-no-special-elig-numpmtct.PJNZ")
demp <- prepare_leapfrog_demp(pjnz_child)
proj <- prepare_leapfrog_projp(pjnz_child)
proj <- prepare_hc_leapfrog_projp(pjnz_child, proj)


demp$netmigr <- leapfrog:::read_netmigr(pjnz_child, adjust_u5mig = FALSE)
demp$netmigr_adj <- leapfrog:::adjust_spectrum_netmigr(demp$netmigr)

dpfile <- grep(".DP$", utils::unzip(pjnz_child, list=TRUE)$Name, value=TRUE)
dp <- utils::read.csv(unz(pjnz_child, dpfile), as.is=TRUE)
dpsub <- function(dp, tag, rows, cols, tagcol=1){
  dp[which(dp[,tagcol]==tag)+rows, cols]
}
yr_start <- as.integer(dpsub(dp,"<FirstYear MV2>",2,4))
yr_end <- as.integer(dpsub(dp, "<FinalYear MV2>",2,4))
proj.years <- yr_start:yr_end
timedat.idx <- 4+1:length(proj.years)-1

pop1 = paste0(getwd(), '/', gsub(x = pjnz_child, pattern = '.PJNZ', replacement = '_pop1.xlsx'))
#pop1 = paste0( gsub(x = pjnz_child, pattern = '.PJNZ', replacement = '_pop1.xlsx'))

spectrum_output <- function(file = "../testdata/spectrum/v6.13/bwa_aim-adult-child-input-art-elig_spectrum-v6.13_2022-02-12_pop1.xlsx", ages = 0:14, country = 'Botswana', years_in = 1970:2030){
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself
  df <- file
  # if(grepl(pattern = 'testdata', file)){
  #   df <- test_path(df)
  # }
  df <- eppasm::read_pop1(df, country, years = years_in)
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
  # setnames(df, 'cd4', 'cd4_cat')
  df <- df %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission)

  df_on_treatment <- df[grepl('ART', df$transmission),] %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission) %>% dplyr::group_by(sex, age, cd4_cat, transmission, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique() %>% dplyr::filter(age %in% ages)
  df_off_treatment <- df[!grepl('ART', df$transmission),]%>% dplyr::filter(age %in% ages)
  df_total <- df %>% dplyr::select(sex, age, cd4_cat, year, pop) %>% dplyr::group_by(sex, age, cd4_cat, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique()%>% dplyr::filter(age %in% ages)

  return(list(on_treatment = df_on_treatment, off_treatment = df_off_treatment, total = df_total))

}
df <- spectrum_output(pop1, ages = 0:14, 'country', years_in = 1970:2030)
x = data.table(df$total)

tag.x ="<AIDSDeathsNoARTSingleAge MV>"
start.id = 20898
end.id = 21148
aids_deathsnoart <- array(as.numeric(unlist(dpsub(dp, tag.x, 3:(end.id - start.id - 2), timedat.idx))), dim = c(length(3:(end.id - start.id - 2)),length(timedat.idx)))
m = aids_deathsnoart[84:98,]
f = aids_deathsnoart[166:180,]
aids_deathsnoart <- array(0, dim = c(15,2,61))
aids_deathsnoart[,1,] <- m
aids_deathsnoart[,2,] <- f

tag.x ="<AIDSDeathsARTSingleAge MV>"
start.id = 20608
end.id = 20858
aids_deathsart <- array(as.numeric(unlist(dpsub(dp, tag.x, 3:(end.id - start.id - 2), timedat.idx))), dim = c(length(3:(end.id - start.id - 2)),length(timedat.idx)))
m = aids_deathsart[84:98,]
f = aids_deathsart[166:180,]
aids_deathsart <- array(0, dim = c(15,2,61))
aids_deathsart[,1,] <- m
aids_deathsart[,2,] <- f

saveRDS(list(proj = proj, demp = demp, dp = dp, timedat.idx = timedat.idx, pjnz = pjnz_child,
             pop1_outputs = x, on_treatment = df$on_treatment, off_trt = df$off_treatment,
             deaths_noart = aids_deathsnoart,
             deaths_art = aids_deathsart),
        testthat::test_path("testdata/child_parms.rds"))
