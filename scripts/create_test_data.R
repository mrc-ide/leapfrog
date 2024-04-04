#!/usr/bin/env Rscript
library(leapfrog)
library(data.table)
library(dplyr)

setwd('C:/Users/mwalters/frogger/')

## Create demographic and projection parameters for adults
pjnz1 <- testthat::test_path("testdata/bwa_aim-no-special-elig.PJNZ")

demp <- prepare_leapfrog_demp(pjnz1)
saveRDS(demp, testthat::test_path("testdata/demographic_projection_object_adult.rds"))
proj <- prepare_leapfrog_projp(pjnz1)
saveRDS(proj, testthat::test_path("testdata/projection_parameters_adult.rds"))

## Used as reference data (Run from leapfrog/master)
# lmod <- leapfrogR(demp, proj)
# saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit.rds"))
#
# lmod <- leapfrogR(demp, proj, hiv_strat = "coarse")
# saveRDS(lmod, testthat::test_path("testdata/leapfrog_fit_coarse.rds"))
#
# mod <- leapfrogR(demp, proj, hiv_steps_per_year = 0L)
# saveRDS(lmod, testthat::test_path("testdata/fit_demography.rds"))


#Create paeds parameters (Run from leapfrog/uncertainrt_analysis_working)
demp <- prepare_leapfrog_demp(pjnz1)
proj <- prepare_leapfrog_projp(pjnz1)
proj <- leapfrog:::prepare_hc_leapfrog_projp(pjnz1, proj)
lmod <- leapfrogR(demp, proj)

demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

setwd('C:/Users/mwalters/frogger/tests/testthat/testdata/')
saveRDS(proj, "projection_parameters_child.rds")
saveRDS(demp, "demographic_projection_object_child.rds")



