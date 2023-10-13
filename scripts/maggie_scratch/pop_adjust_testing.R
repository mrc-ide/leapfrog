rm(list = ls())
library(tidyr)
library(tidyverse)
library(dplyr)
library(data.table)
setwd('C:/Users/mwalters/leapfrog/')
devtools::load_all()
pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")
demp <- readRDS(("C:/Users/mwalters/frogger/tests/testthat/testdata/demographic_projection_object_child.rds"))
parameters <- readRDS(("C:/Users/mwalters/frogger/tests/testthat/testdata/projection_parameters_child.rds"))


parameters$ctx_effect[] <- 0
parameters$ctx_val[] <- 0
save = parameters$pmtct
parameters$pmtct <- parameters$pmtct[,,2]
parameters$pmtct_input_isperc[] <- as.integer(1)
parameters$artpaeds_isperc[] <- TRUE
parameters$paed_art_val[] <- 0
parameters$pmtct_mtct[] <- 0
parameters$mtct[] <- 0
demp$targetpop <- demp$basepop
out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE, pop_adjust = T)


basepop <- data.table(melt(demp$basepop))
basepop$Var2 <- as.integer(basepop$Var2)
totalpop <- data.table(melt(out$p_total_pop))
totalpop <- totalpop[,.(Var1 = Var1 - 1, Var2, Var3 = Var3 + 1969, value)]
pop <- merge(basepop, totalpop, by = c('Var1', 'Var2', 'Var3'))
pop[,diff := value.x - value.y]


fr <- data.table(melt(out$h_hiv_adult))
fr <- fr[,.(cd4 = Var1, age = Var2 + 14, sex = ifelse(Var3 == 1, 'Male', 'Female'), year = Var4 + 1969, fr = value)]
spec = spectrum_output(file = 'C:/Users/mwalters/leapfrog/tests/testdata/spectrum/v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18_pop1.xlsx', ages = 15:49, country = 'bwa')

spec <- data.table(spec$total)

cd4_cat = data.table(cd4 = 1:7, cd4_cat = c('gte500', '350-500', '250-349', '200-249', '100-199','50-99', 'lte50'))
fr <- merge(fr, cd4_cat, by = 'cd4')
fr[,cd4 := NULL]

dt <- merge(spec, fr, by = c('age', 'sex', 'year', 'cd4_cat'))
