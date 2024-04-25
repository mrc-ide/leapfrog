#!/usr/bin/env Rscript

## This script creates test data required for the frogger tests.
## We read some input data and prepare a set of demographic projection
## and HIV parameters for both the adult and the child model.
## We also run leapfrog and save out the result for use in reference tests

# nolint start
library(leapfrog)
library(data.table)
library(dplyr)

## Create demographic and projection parameters for adults
pjnz1 <- testthat::test_path("testdata/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")

demp <- prepare_leapfrog_demp(pjnz1)
saveRDS(demp, testthat::test_path("testdata/demographic_projection_object_adult.rds"))
proj <- prepare_leapfrog_projp(pjnz1)
saveRDS(proj, testthat::test_path("testdata/projection_parameters_adult.rds"))

## Used as reference data
lmod <- leapfrogR(demp, proj)
saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit.rds"))

lmod <- leapfrogR(demp, proj, hiv_strat = "coarse")
saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit_coarse.rds"))

mod <- leapfrogR(demp, proj, hiv_steps_per_year = 0L)
saveRDS(lmod, testthat::test_path("testdata/fit_demography.rds"))


## TODO: this should be built from the PJNZ source data
# pjnz1 <- testthat::test_path("testdata/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")

hivp <- readRDS(testthat::test_path("testdata/projection_parameters_child.rds"))

## We need 61 entries for 1970 to 2030 inclusive. If the input data
## has fewer entries than this, then add an additional value
if (length(hivp$artpaeds_isperc) < 61) {
  hivp$artpaeds_isperc <- c(hivp$artpaeds_isperc, FALSE)
}
if (length(hivp$paed_art_elig_age) < 61) {
  hivp$paed_art_elig_age <- c(hivp$paed_art_elig_age, 2)
  hivp$paed_art_elig_age <- as.integer(hivp$paed_art_elig_age)
}
if (ncol(hivp$paed_art_elig_cd4) < 61) {
  hivp$paed_art_elig_cd4 <- cbind(hivp$paed_art_elig_cd4, hivp$paed_art_elig_cd4[,ncol(hivp$paed_art_elig_cd4)])
}

hivp$laf = 1

saveRDS(hivp, testthat::test_path("testdata/projection_parameters_child.rds"))
# nolint end
