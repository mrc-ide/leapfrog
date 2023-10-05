rm(list = ls())
library(tidyr)
library(tidyverse)
library(dplyr)
setwd('C:/Users/mwalters/leapfrog/')
devtools::load_all()
pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")
demp <- readRDS(("C:/Users/mwalters/frogger/tests/testthat/testdata/demographic_projection_object_child.rds"))
parameters <- readRDS(("C:/Users/mwalters/frogger/tests/testthat/testdata/projection_parameters_child.rds"))


parameters$ctx_effect[] <- 0
parameters$ctx_val[] <- 0
# parameters$mtct[] <- 0
#parameters$pmtct[] <- 0
save = parameters$pmtct
parameters$pmtct <- parameters$pmtct[,,2]
# parameters$pmtct_dropout[] <- 0
parameters$pmtct_input_isperc[] <- as.integer(1)
# parameters$bf_duration_art[] <- 1
# parameters$bf_duration_no_art[] <- 1
parameters$incidinput[which(1970:2030 == 1991):length(1970:2030)] <- 0
parameters$paed_incid_input[which(1970:2030 == 2000)] <- 0
# parameters$bf_duration[] <- 1
parameters$paed_art_elig_age <- as.integer(parameters$paed_art_elig_age)
out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)

parameters$pmtct <- save
parameters$pmtct[,,2] <- parameters$pmtct[,,2]/ 100
parameters$paed_art_elig_age <- as.numeric(parameters$paed_art_elig_age)
parameters$pmtct_input_isperc <- as.logical(parameters$pmtct_input_isperc[])
pmtct_mtct <- array(unlist(list(cbind(parameters$mtct[,1], parameters$pmtct_mtct[,,1]), cbind(parameters$mtct[,2], parameters$pmtct_mtct[,,2]))), dim = c(7,8,2))
parameters$pmtct_mtct <- pmtct_mtct
lmod <- leapfrogR(demp, parameters)
save <- (lmod)
library(data.table)

# if(time_step == 26){
#   std::cout << intermediate.children.prop_wlhiv_gte350;
# }

# #
# # ##Check number of women with HIV
fr <- out$h_hiv_adult[,1:35,2,]
fr <- apply(fr, MARGIN = length(dim(fr)), FUN = sum)
lf <- lmod$hivstrat_adult[,1:35,2,]
lf <- apply(lf, MARGIN = length(dim(lf)), FUN = sum)
# plot(x = 1970:2030, y = fr) + lines(x = 1970:2030, y = lf)
dt <- data.table(year = 1970:2030, lfrog = lf, frogger = fr)
dt[,diff := lfrog - frogger]
dt

fr_art <- out$h_art_adult[,,1,2,]
lf_art <- lmod$artstrat_adult[,,1,2,]
fr_art <- apply(fr_art, MARGIN = 3, FUN = sum)
lf_art <- apply(lf_art, MARGIN = 3, FUN = sum)
dt_art <- data.table(year = 1970:2030, lfrog = lf_art, frogger = fr_art)
dt_art[,diff := lfrog - frogger]
dt_art


###Check number of births to WLHIV
fr.x = as.vector(out$hiv_births)
lfrog.x = as.vector(lmod$hiv_births)
births = data.table(year = 1970:2030, fr = fr.x, lfrog = lfrog.x)
births[,diff := fr - lfrog]
births

##NO ART
{
lmod = save
l <- data.table(melt(lmod$hivstrat_paeds[1:6, , 6:15, , ]))
o <- data.table(melt(out$hc2_hiv_pop))

setnames(o, 'value', 'value_o')
dt <- merge(l, o, by = c('Var1', 'Var2', 'Var3', 'Var4', 'Var5'))
dt[,diff := value - value_o]
dt2 <- dt[,.(cd4 = Var1, cat = Var2, age = Var3 + 4, sex = Var4, year = Var5 + 1969, value, value_o, diff)]

l <- data.table(melt(lmod$hivstrat_paeds[1:7, , 1:5, , ]))
o <- data.table(melt(out$hc1_hiv_pop))

setnames(o, 'value', 'value_o')
dt <- merge(l, o, by = c('Var1', 'Var2', 'Var3', 'Var4', 'Var5'))
dt[,diff := value - value_o]
dt1 <- dt[,.(cd4 = Var1, cat = Var2, age = Var3 - 1, sex = Var4, year = Var5 + 1969, value, value_o, diff)]

dt <- rbind(dt1, dt2)


noart <- dt[,type := 'noart']
}

##ART
{
  save = lmod
  l <- data.table(melt(lmod$artstrat_paeds[, 1:6, 6:15, , ]))
  o <- data.table(melt(out$hc2_art_pop))


  setnames(o, 'value', 'value_o')
  dt <- merge(l, o, by = c('Var1', 'Var2', 'Var3', 'Var4', 'Var5'))
  dt[,diff := value - value_o]
  dt2 <- dt[,.(trt_time = Var1, cd4  = Var2, age = Var3 + 4, sex = Var4, year = Var5 + 1969, value, value_o, diff)]

  l <- data.table(melt(lmod$artstrat_paeds[, , 1:5, , ]))
  o <- data.table(melt(out$hc1_art_pop))

  setnames(o, 'value', 'value_o')
  dt <- merge(l, o, by = c('Var1', 'Var2', 'Var3', 'Var4', 'Var5'))
  dt[,diff := value - value_o]
  dt1 <- dt[,.(trt_time = Var1, cd4 = Var2, age = Var3 - 1, sex = Var4, year = Var5 + 1969, value, value_o, diff)]

  dt <- rbind(dt1, dt2)

  art = dt[,type := 'art']
  }


dt <- rbind(art, noart, fill = T)
x = dt[value_o > 1e10]
x[order(year)]


##Focus on ones that aren't aligning
prob = dt[abs(diff) > 1e-3]




##deaths
{
deaths <- data.table(melt((out$hc2_noart_aids_deaths)))
  deaths <- deaths[,.(cd4 = Var1, tt = Var2, age = Var3 + 4, sex = Var4, year = Var5 + 1969, value)]

  deaths_lmod <- data.table(melt(lmod$deaths_paeds))
  deaths_lmod <- deaths_lmod[,.(cd4 = Var1, tt = Var2, age = Var3 - 1, sex = Var4, year = Var5 + 1969, value)]


}



