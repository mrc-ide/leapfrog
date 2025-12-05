#!/usr/bin/env Rscript

## This script creates test data required for the leapfrog tests.
## We read some input data and prepare a set of demographic projection
## and HIV parameters for both the adult and the child model.
## We also run leapfrog and save out the result for use in reference tests
setwd("leapfrogr")
devtools::load_all()

# nolint start
library(dplyr)

## Create demographic and projection parameters for adults
pjnz_adult <- file.path(here::here(), "inst", "pjnz", "bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")

parameters <- process_pjnz(pjnz_adult)
save_parameters(parameters, testthat::test_path("testdata/adult_parms_full.h5"))

parameters_coarse <- process_pjnz(pjnz_adult, use_coarse_age_groups = TRUE)
save_parameters(parameters_coarse, testthat::test_path("testdata/adult_parms_coarse.h5"))

# We use France for testing Spectrum model variant as it has non-zero Non-AIDS excess mortality
# inputs. Which we need for testing modelled non-AIDS excess mortality.
# This was created by creating a new default projection in Spectrum
pjnz_france <- file.path(here::here(), "inst", "pjnz", "france_default.PJNZ")
france_parameters <- process_pjnz(pjnz_france)
save_parameters(france_parameters, testthat::test_path("testdata/spectrum_params.h5"))

#Create paeds parameters
pjnz_child <- file.path(here::here(), "inst", "pjnz", "bwa_aim-no-special-elig-numpmtct.PJNZ")

parameters <- process_pjnz(pjnz_child)
parameters_coarse <- process_pjnz(pjnz_child, use_coarse_age_groups = TRUE)

pop1 <- gsub(pattern = '.PJNZ', replacement = '_pop1.xlsx', x = pjnz_child)

spectrum_output <- function(file, ages = 0:14, country = 'Botswana', years_in = 1970:2030){
  ##pull out stratified population from the .xlsx file, This function doesn't take out the paediatric output, so going to just compare to the Spectrum software itself
  df <- file
  # if(grepl(pattern = 'testdata', file)){
  #   df <- test_path(df)
  # }
  df <- eppasm::read_pop1(df, country, years = years_in)
    df_paed <- df %>% dplyr::filter(age < 5) %>%
      dplyr::right_join(y = data.frame(cd4 = 1:8, cd4_cat = c('neg', 'gte30', '26-30', '21-25', '16-20', '11-15', '5-10', 'lte5'))) %>%
      dplyr::right_join(y = data.frame(artdur = 2:8, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg')

    df_adol <- df %>% dplyr::filter(age > 4 & age < 15) %>%
      dplyr::right_join(y = data.frame(cd4 = 3:8, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200'))) %>%
      dplyr::right_join(y = data.frame(artdur = 2:8, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg')

    df_paed <- rbind(df_paed, df_adol)

    df <- df %>% dplyr::filter(age %in%  c(15:max(ages))) %>%
      dplyr::right_join(y = data.frame(cd4 = 2:8, cd4_cat = c('gte500', '350-500', '250-349', '200-249', '100-199','50-99', 'lte50'))) %>%
      dplyr::right_join(y = data.frame(artdur = 2:8, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'ARTlte5mo', 'ART6to12mo', 'ARTgte12mo'))) %>%
      dplyr::filter(cd4_cat != 'neg')

    df <- rbind(df, df_paed)
  # setnames(df, 'cd4', 'cd4_cat')
  df <- df %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission)

  df_on_treatment <- df[grepl('ART', df$transmission),] %>% dplyr::select(sex, age, cd4_cat, year, pop, transmission) %>% dplyr::group_by(sex, age, cd4_cat, transmission, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique() %>% dplyr::filter(age %in% ages)
  df_off_treatment <- df[!grepl('ART', df$transmission),]%>% dplyr::filter(age %in% ages)
  df_total <- df %>% dplyr::select(sex, age, cd4_cat, year, pop) %>% dplyr::group_by(sex, age, cd4_cat, year) %>% dplyr::mutate(pop = sum(pop)) %>% unique()%>% dplyr::filter(age %in% ages)

  return(list(on_treatment = df_on_treatment, off_treatment = df_off_treatment, total = df_total))

}
df <- spectrum_output(pop1, ages = 0:80, 'country', years_in = 1970:2030)
x <- df$total

dp <- read_dp(pjnz_child)
dat <- parse_dp(dp)

aids_deathsnoart <- array(0, dim = c(15,2,61), dimnames = list(age = 0:14, sex = c('male', 'female'), years = 1970:2030))
aids_deathsnoart[,'male',] <- dat$data$aids_deaths_no_art_single_age$data[as.character(0:14),"male",]
aids_deathsnoart[,'female',] <- dat$data$aids_deaths_no_art_single_age$data[as.character(0:14),"female",]

aids_deathsart <- array(0, dim = c(15,2,61), dimnames = list(age = 0:14, sex = c('male', 'female'), years = 1970:2030))
aids_deathsart[,'male',] <- dat$data$aids_deaths_art_single_age$data[as.character(0:14),"male",]
aids_deathsart[,'female',] <- dat$data$aids_deaths_art_single_age$data[as.character(0:14),"female",]

spec_ctx_need <- dat$data$child_art_calc$data["both", "Children needing cotrim (0-14): ",]

out <- list(dp = dp,
            pjnz = pjnz_child,
            pop1 = x,
            ontrt = df$on_treatment,
            offtrt = df$off_treatment,
            deaths_noart = aids_deathsnoart,
            deaths_art = aids_deathsart,
            ctx_need = as.numeric(unlist(spec_ctx_need)))

save_parameters(parameters, testthat::test_path("testdata/child_parms_full.h5"))
save_parameters(parameters_coarse, testthat::test_path("testdata/child_parms_coarse.h5"))

# need this for child model tests
saveRDS(out, testthat::test_path("testdata/child_test_utils.rds"))
