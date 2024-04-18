#!/usr/bin/env Rscript
library(leapfrog)
library(data.table)
library(dplyr)
devtools::load_all('C:/Users/mwalters/leapfrog')
setwd('C:/Users/mwalters/frogger/')

## Create demographic and projection parameters for adults
pjnz1 <- testthat::test_path("testdata/bwa_aim-no-special-elig-numpmtct.PJNZ")
# pjnz1 <- testthat::test_path("testdata/bwa_aim-no-special-elig_ctx.PJNZ")
#pjnz1 <- testthat::test_path("testdata/bwa_aim-no-special-elig.PJNZ")
#pjnz1 = 'C:/Users/mwalters/Desktop/NW_TEST_MTCT_BF_PERI_pmtct.PJNZ'

demp <- prepare_leapfrog_demp(pjnz1)
saveRDS(demp, testthat::test_path("testdata/demographic_projection_object_adult.rds"))
proj <- prepare_leapfrog_projp(pjnz1)
saveRDS(proj, testthat::test_path("testdata/projection_parameters_adult.rds"))

# Used as reference data (Run from leapfrog/master)
lmod <- leapfrogR(demp, proj)
saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit.rds"))

lmod <- leapfrogR(demp, proj, hiv_strat = "coarse")
saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit_coarse.rds"))

mod <- leapfrogR(demp, proj, hiv_steps_per_year = 0L)
saveRDS(lmod, testthat::test_path("testdata/fit_demography.rds"))


#Create paeds parameters (Run from leapfrog/uncertainrt_analysis_working)
demp <- prepare_leapfrog_demp(pjnz1)
proj <- prepare_leapfrog_projp(pjnz1)
proj <- leapfrog:::prepare_hc_leapfrog_projp(pjnz1, proj)


demp$netmigr <- leapfrog:::read_netmigr(pjnz1, adjust_u5mig = FALSE)
demp$netmigr_adj <- leapfrog:::adjust_spectrum_netmigr(demp$netmigr)

dpfile <- grep(".DP$", utils::unzip(pjnz1, list=TRUE)$Name, value=TRUE)
dp <- utils::read.csv(unz(pjnz1, dpfile), as.is=TRUE)
dpsub <- function(tag, rows, cols, tagcol=1){
  dp[which(dp[,tagcol]==tag)+rows, cols]
}
yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
proj.years <- yr_start:yr_end
timedat.idx <- 4+1:length(proj.years)-1

pop1 = paste0(getwd(), '/', gsub(x = pjnz1, pattern = '.PJNZ', replacement = '_pop1.xlsx'))
#pop1 = paste0( gsub(x = pjnz1, pattern = '.PJNZ', replacement = '_pop1.xlsx'))

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

saveRDS(list(proj = proj, demp = demp, dp = dp, timedat.idx = timedat.idx, pjnz = pjnz1,
             pop1_outputs = x, on_treatment = df$on_treatment, off_trt = df$off_treatment),
        "C:/Users/mwalters/frogger/tests/testthat/testdata/child_parms.rds")



