#!/usr/bin/env Rscript
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


# TODO: Add details about child model input data
#
# Test data coming from the same file but on the 'clean' branch of leapfrog.
#
# ```R
# pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")
#
# demp <- prepare_leapfrog_demp(pjnz1)
# hivp <- prepare_leapfrog_projp(pjnz1)
# hivp = prepare_hc_leapfrog_projp(pjnz1, hivp)
#
# demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
# demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
# hivp$incidinput[which(1970:2030 == 2000):length(parameters$incidinput)] <- 0
#
# hivp$paed_cd4_dist <- c(1,rep(0,6))
#
# ##HIV starts in 2000
# hivp$paed_incid_input[] <- 0
# hivp$paed_incid_input[which(1970:2030 ==2000)] <- 100
#
# ##Only transmission is coming from nosocomial infections
# hivp$pmtct_mtct[] <- 0
# hivp$art_mtct[] <- 0
#
# ##Add in treatment
# hivp$paed_art_val[which(1970:2030 == 2002)] <- 50
# hivp$artpaeds_isperc[] <- FALSE
#
# ##Change things to length 61
# hivp$artpaeds_isperc <- c(hivp$artpaeds_isperc, FALSE)
# hivp$paed_art_elig_age <- c(hivp$paed_art_elig_age, 2)
# hivp$paed_art_elig_cd4 <- cbind(hivp$paed_art_elig_cd4, hivp$paed_art_elig_cd4[,ncol(hivp$paed_art_elig_cd4)])
#
# hivp$laf = 1
#
#
# setwd('C:/Users/mwalters/frogger/tests/testthat/testdata/')
# saveRDS(hivp, "projection_parameters_child.rds")
# saveRDS(demp, "demographic_projection_object_child.rds")
#
# ```

# nolint end
