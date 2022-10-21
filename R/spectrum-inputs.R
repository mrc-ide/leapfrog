
read_sx <- function(pjnz, use_ep5=FALSE){

  if(use_ep5) {
    dpfile <- grep(".ep5$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  } else {
    dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  }

  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  ## projection parameters
  yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1

  ## mx
  Sx <- dpsub("<SurvRate MV2>", 3+c(0:81, 82+0:81), timedat.idx)
  Sx <- array(as.numeric(unlist(Sx)), c(82, 2, length(proj.years)))
  dimnames(Sx) <- list(age=c(0:80, "80+"), sex=c("male", "female"), year=proj.years)

  Sx
}

read_netmigr <- function(pjnz, use_ep5=FALSE, adjust_u5mig = TRUE, sx = NULL){

  if(use_ep5) {
    dpfile <- grep(".ep5$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  } else {
    dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  }

  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  ## projection parameters
  yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1

  ## netmig
  totnetmig <- sapply(dpsub("<MigrRate MV2>", c(4, 6), timedat.idx), as.numeric)

  netmigagedist <- sapply(dpsub("<MigrAgeDist MV2>", 2 + 1:34, timedat.idx), as.numeric)/100
  netmigagedist <- array(c(netmigagedist), c(17, 2, length(proj.years)))
  netmigr5 <- sweep(netmigagedist, 2:3, totnetmig, "*")

  netmigr <- array(dim = c(81, 2, length(proj.years)),
                   dimnames = list(age = 0:80, sex = c("male", "female"), year = proj.years))
  netmigr[1:80, , ] <- apply(netmigr5[1:16,,], 2:3, beers::beers_sub_ordinary)
  netmigr[81, , ] <- netmigr5[17, , ]

  if (adjust_u5mig) {

    if (is.null(sx)) {
      sx <- read_sx(pjnz)
    }
    
    u5prop <- array(dim = c(5, 2))
    u5prop[1, ] <- sx[1, , 1] * 2
    u5prop[2, ] <- sx[2, , 1] * u5prop[1, ]
    u5prop[3, ] <- sx[3, , 1] * u5prop[2, ]
    u5prop[4, ] <- sx[4, , 1] * u5prop[3, ]
    u5prop[5, ] <- sx[5, , 1] * u5prop[4, ]
    
    u5prop <- sweep(u5prop, 2, colSums(u5prop), "/")
    
    netmigr[1:5 , 1, ] <- u5prop[ , 1, drop = FALSE] %*% netmigr5[1, 1, ]
    netmigr[1:5 , 2, ] <- u5prop[ , 2, drop = FALSE] %*% netmigr5[1, 2, ]
  }
  
  netmigr
}

adjust_spectrum_netmigr <- function(netmigr) {

  ## Spectrum adjusts net-migration to occur half in
  ## current age group and half in next age group
  
  netmigr_adj <- netmigr
  netmigr_adj[-1,,] <- (netmigr[-1,,] + netmigr[-81,,])/2
  netmigr_adj[1,,] <- netmigr[1,,]/2
  netmigr_adj[81,,] <- netmigr_adj[81,,] + netmigr[81,,]/2

  netmigr_adj
}


#' Prepare demographic inputs from Spectrum PJNZ
#'
#' @param pjnz path to PJNZ file
#'
#' @return list of demographic input parameters
#'
#' @examples
#' pjnz <- system.file("pjnz/bwa2021_v6.13.pjnz", package = "leapfrog")
#' demp <- prepare_leapfrog_demp(pjnz)
#'
#' @export
prepare_leapfrog_demp <- function(pjnz) {

  demp <- eppasm::read_specdp_demog_param(pjnz)
  demp$Sx <- read_sx(pjnz)
  demp$netmigr <- read_netmigr(pjnz, sx = demp$Sx)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

  demp$births_sex_prop <- rbind(male = demp$srb, female = 100) / (demp$srb + 100)

  ## normalise ASFR distribution
  demp$asfr <- sweep(demp$asfr, 2, demp$tfr / colSums(demp$asfr), "*")
  

  demp
}

#' Prepare adult HIV projection parameters from Spectrum PJNZ
#'
#' @param pjnz path to PJNZ file
#' @param hiv_steps_per_year number of Euler integration steps per year; default 10
#' @param hTS number of HIV treatment stages; default 3 (0-5 months,
#'   6-11 months, 12+ months)
#'
#' @return list of HIV projection parameters
#'
#' @examples
#' pjnz <- system.file("pjnz/bwa2021_v6.13.pjnz", package = "leapfrog")
#' projp <- prepare_leapfrog_projp(pjnz)
#'
#' @export
prepare_leapfrog_projp <- function(pjnz, hiv_steps_per_year = 10L, hTS = 3) {
 
  projp <- read_hivproj_param(pjnz)
    
  ## Hard coded to expand age groups 15-24, 25-34, 35-44, 45+ to
  ## single-year ages 15:80.
  ## Requires extension for coarse HIV age group stratification
  idx_expand_full <- rep(1:4, times = c(10, 10, 10, 36))
  idx_expand_coarse <- rep(1:4, times = c(3, 2, 2, 2))

  v <- list()
  v$incidinput <- eppasm::read_incid_input(pjnz)
  v$incidpopage <- attr(v$incidinput, "incidpopage")
  v$incrr_sex <- projp$incrr_sex
  
  ## HIV effects on fertilty
  ## mkw: should this be normalized as well?
  x <- apply(projp$fert_rat, 2, rep, each = 5)
  x <- rbind(array(data = 0, dim = c(15,61)), x, array(data = 0, dim = c(31,61)))
  rownames(x) <- 0:80
  v$fert_rat <- x
  
  ## paed input
  v$paed_incid_input <- projp$nosocom_infections_04
  ## Hardcoded, this is putting all individuals in the highest cd4 category bc i think thats how the nosocomial infections work
 ## v$paed_cd4_dist <- c(0.6, 0.12, 0.1, 0.09, 0.05, 0.03, 0.01)
  v$paed_cd4_dist <- c(1, 0, 0, 0, 0, 0, 0)
  paed_cd4_prog <- array(data = 0, dim = c(6, 5, 2), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10'), age = 0:4, sex = c('male', 'female')) )
  paed_cd4_prog[,1:5,1:2] <- c(0.14, 0.37, 0.3, 0.35, 0.4, 0.4)
  
  adol_cd4_prog <- array(data = 0, dim = c(5, 10, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349'), age = 5:14 , sex = c('male', 'female')) )
  adol_cd4_prog[,1:10,1:2] <- c(0.3028, 0.3028, 0.2575, 0.2122, 0.1669)
  
  v$paed_cd4_prog <- paed_cd4_prog
  v$adol_cd4_prog <- adol_cd4_prog
  
  paed_cd4_mort <- array(data = 0, dim = c(7, 4, 5, 2), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                        transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                        age = c(0:4), sex = c('male', 'female')))
  ## 0-2
  paed_cd4_mort[,1,1:3,] <- c(0.25801, 0.31513, 0.38490, 0.47011, 0.57418, 0.79129, 0.85655)
  paed_cd4_mort[,2,1:3,] <- c(0.15385, 0.18791, 0.22951, 0.28031, 0.34237, 0.41817, 0.51074)
  paed_cd4_mort[,3,1:3,] <- rep(0.08743, times = 7)
  paed_cd4_mort[,4,1:3,] <- rep(0.02450, times = 7)
  
  ## 3-4
  paed_cd4_mort[,1,4:5,] <- c(0.06470, 0.07902, 0.09652, 0.11789, 0.14398, 0.14586, 0.21479)
  paed_cd4_mort[,2,4:5,] <- rep(0.04622, times = 7)
  paed_cd4_mort[,3,4:5,] <- c(0.03383, 0.03485, 0.03609, 0.03761, 0.03946, 0.04172, 0.04449)
  paed_cd4_mort[,4,4:5,] <- c(0.02144, 0.02347, 0.02596, 0.02899, 0.03270, 0.03723, 0.04276)
  
  
  paed_art_mort <- array(data = 0, dim = c(7, 3, 5, 2), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                        transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                        age = c(0:4), sex = c('male', 'female')))
  ## 0-6 mo on treatment
  paed_art_mort[,1,1,] <- c(0.133460, 0.142490, 0.154150, 0.170110, 0.194100, 0.237140, 0.364790)
  paed_art_mort[,1,2:3,] <- c(0.061010, 0.065140, 0.070470, 0.077770, 0.088730, 0.108400, 0.166760)
  paed_art_mort[,1,4:5,] <- c(0.028720, 0.030660, 0.033170, 0.036600, 0.041760, 0.051020, 0.078490)

  ## 7-12 mo on treatment
  paed_art_mort[,2,1,] <- c(0.048660, 0.051950, 0.056200, 0.062020, 0.070770, 0.086460, 0.133000)
  paed_art_mort[,2,2:3,] <- c(0.022240, 0.023750, 0.025690, 0.028350, 0.032350, 0.039520, 0.060800)
  paed_art_mort[,2,4:5,] <- c(0.010470, 0.011180, 0.012090, 0.013340, 0.015230, 0.018600, 0.028620)
  
  ## 12+ mo on treatment
  paed_art_mort[,3,1,] <- rep(0, 7)
  paed_art_mort[,3,2:3,] <- c(0.018610, 0.019340, 0.020260, 0.021480, 0.023220, 0.026140, 0.033720)
  paed_art_mort[,3,4:5,] <- c(0.008760, 0.009100, 0.009540, 0.010110, 0.010930, 0.012300, 0.015870)
  
  adol_cd4_mort <- array(data = 0, dim = c(6, 4, 10, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                        transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                        age =5:14, sex = c('male', 'female')))
  ## 5 - 14
  adol_cd4_mort[,1,,] <- c(0.02617, 0.03197, 0.03905, 0.04769, 0.05825, 0.07114)
  adol_cd4_mort[,2,,] <- c(0.02572, 0.03142, 0.03837, 0.04687, 0.05724, 0.06991)
  adol_cd4_mort[,3,,] <- c(0.02409, 0.02942, 0.03593, 0.04388, 0.05360, 0.06547)
  adol_cd4_mort[,4,,] <- c(0.02245, 0.02742, 0.03349, 0.04090, 0.04996, 0.06102)
  
  adol_art_mort <- array(data = 0, dim = c(6, 3, 10, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                        transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                        age = c(5:14), sex = c('male', 'female')))
  ## 0-6 mo on treatment
  adol_art_mort[,1,1:5,] <- c(0.013800, 0.015860, 0.019110, 0.023650, 0.030090, 0.052640)
  adol_art_mort[,1,6:10,] <- c(0.015140, 0.017400, 0.020950, 0.025930, 0.032990, 0.057730)

  ## 7-12 mo on treatment
  adol_art_mort[,2,1:5,] <- c(0.005180, 0.005960, 0.007170, 0.008880, 0.011300, 0.019770)
  adol_art_mort[,2,6:10,] <- c(0.005680, 0.006530, 0.007870, 0.009740, 0.012390, 0.021680)

  ## 12+ mo on treatment
  adol_art_mort[,3,1:5,] <- c(0.005240, 0.005620, 0.006170, 0.006880, 0.007770, 0.010320)
  adol_art_mort[,3,6:10,] <- c(0.005740, 0.006160, 0.006770, 0.007540, 0.008520, 0.011320)
  
  v$paed_cd4_mort <- paed_cd4_mort
  v$adol_cd4_mort <- adol_cd4_mort
  v$paed_art_mort <- paed_art_mort
  v$adol_art_mort <- adol_art_mort
  
  ## pull in cotrim coverage numbers
  ctx_val <- c(rep(0, 25), rep(0, 5), 
               rep(50, 3),
               52.222220,	54.444440,	56.666670,	58.888890,	61.111110,	65.000000,	68.888890,	72.777780,	76.666670,
               80.555550,	84.444440,	88.333330,	92.222220,	96.111110,	86.690000,	100.000000,	100.000000,	100.000000,
               100.000000,	100.000000,	100.000000,	100.000000,	100.000000,	100.000000,	100.000000,	100.000000,	100.000000,	100.000000)
  ctx_pct <- T
  if(ctx_pct){
    ctx_val <- ctx_val / 100
  }
 ## ctx_effect_notrt <- c(rep(0.33, 5), rep(0,5))
  ctx_effect <- 0.33
  v$ctx_val <- ctx_val
  v$ctx_effect <- ctx_effect
  
  ## pull in ART coverage numbers
  paed_art_val <- c(rep(0, length(1970:1994)), 0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	13.407190,	27.990470,	44.005150,	58.008330,	69.103150,	79.204180,	75.471110,	72.889590,	80.878660,	89.824530,	94.705560,	86.893400,	88.066360,	72.887590,	59.954710,	59.831570,	50.672620,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990,	62.719990)
  paed_art_pct <- T
  if(paed_art_pct){
    paed_art_val <- paed_art_val / 100
  }
  
  v$paed_art_val <- paed_art_val
  
  ## ART eligibility age, doing this in years rather than months
  v$paed_art_elig_age <- c(rep(0, 37), rep(1, 3), rep(2, 20))
  paed_art_elig_cd4 <- array(data = NA, dim = c(length(0:14), length(1970:2029)), dimnames = list(age = c(0:14), year = c(1970:2029)))
  ## corresponds with less than 25
  paed_art_elig_cd4[1:5,] <- 3
  ## correspond with 200 then 300 
  paed_art_elig_cd4[6:15,1:40] <- 6
  paed_art_elig_cd4[6:15,41:60] <- 5
  v$paed_art_elig_cd4 <- paed_art_elig_cd4
  
  ## HIV positive entrants, right now just doing those without ART
  v$age15hivpop <- projp$age15hivpop

  ## Use Beer's coefficients to distribution IRRs by age/sex
  Amat <- eppasm:::create_beers(17)
  v$incrr_age <- apply(projp$incrr_age, 2:3, function(x)  Amat %*% x)[16:81, , ] ## !! Hard coded
  v$incrr_age[v$incrr_age < 0] <- 0

  v$cd4_initdist_full <- projp$cd4_initdist[ , idx_expand_full, ]
  v$cd4_prog_full <- (1-exp(-projp$cd4_prog[ , idx_expand_full, ] / hiv_steps_per_year)) * hiv_steps_per_year
  v$cd4_mort_full <- projp$cd4_mort[ ,idx_expand_full, ]
  v$art_mort_full <- projp$art_mort[c(1, 2, rep(3, hTS - 2)), , idx_expand_full, ]

  v$cd4_initdist_coarse <- projp$cd4_initdist[ , idx_expand_coarse, ]
  v$cd4_prog_coarse <- (1-exp(-projp$cd4_prog[ , idx_expand_coarse, ] / hiv_steps_per_year)) * hiv_steps_per_year
  v$cd4_mort_coarse <- projp$cd4_mort[ ,idx_expand_coarse, ]
  v$art_mort_coarse <- projp$art_mort[c(1, 2, rep(3, hTS - 2)), , idx_expand_coarse, ]
  
  paed_cd4_transition <- array(0, dim = c(6,7), dimnames = list(cd4_count = c('gte1000', '750-1000', '500-749', '350-499', '200-349', 'lte200'), cd4_pct = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5')))
  paed_cd4_transition[1:2,1] <- c(0.71, 0.29)
  paed_cd4_transition[2:3,2] <- c(0.6, 0.4)
  paed_cd4_transition[3:4,3] <- c(0.83, 0.17)
  paed_cd4_transition[4:5,4] <- c(0.77, 0.23)
  paed_cd4_transition[5:6,5] <- c(0.89, 0.11)
  paed_cd4_transition[6,6:7] <- c(0.01, 0.01)
  v$paed_cd4_transition <- paed_cd4_transition


  v$artmx_timerr <- projp$artmx_timerr[c(1, 2, rep(3, hTS - 2)), ]

  ## ## ART eligibility and numbers on treatment

  v$art15plus_num <- projp$art15plus_num
  v$art15plus_isperc <- projp$art15plus_numperc == 1

  ## convert percentage to proportion
  v$art15plus_num[v$art15plus_isperc] <- v$art15plus_num[v$art15plus_isperc] / 100


  ## eligibility starts in projection year idx
  ## ## !! NOTE: from EPP-ASM; not yet implemented
  ## v$specpop_percelig <- rowSums(with(projp$artelig_specpop[-1,], mapply(function(elig, percent, year) rep(c(0, percent*as.numeric(elig)), c(year - proj_start, proj_end - year + 1)), elig, percent, year)))
  v$artcd4elig_idx <- findInterval(-projp$art15plus_eligthres, -c(999, 500, 350, 250, 200, 100, 50))

  ## Update eligibility threshold from CD4 <200 to <250 to account for additional
  ## proportion eligible with WHO Stage 3/4.
  v$artcd4elig_idx <- replace(v$artcd4elig_idx, v$artcd4elig_idx==5L, 4L)

  v$pw_artelig <- with(projp$artelig_specpop["PW",], rep(c(0, elig), c(year - projp$yr_start, projp$yr_end - year + 1)))  # are pregnant women eligible (0/1)

  ## ## !! NOTE: from EPP-ASM; not yet implemented
  ## ## percentage of those with CD4 <350 who are based on WHO Stage III/IV infection
  ## v$who34percelig <- who34percelig

  v$art_dropout <- projp$art_dropout/100

  proj_years <- as.integer(projp$yr_end - projp$yr_start + 1L)
  v$t_ART_start <- min(c(unlist(apply(v$art15plus_num > 0, 1, which)), proj_years))

  ## New ART patient allocation options
  v$art_alloc_method <- projp$art_alloc_method
  v$art_alloc_mxweight <- projp$art_prop_alloc[1]

  ## Scale mortality among untreated population by ART coverage
  v$scale_cd4_mort <- projp$scale_cd4_mort

  ## State space dimensions
  v$hAG_SPAN_full <- rep(1L, 66L)
  v$hAG_SPAN_coarse <- c(2L, 3L, 5L, 5L, 5L, 5L, 5L, 5L, 31L)
  
  ## Add in pediatric components
  v$fert_rat <- apply(projp$fert_rat, 2, rep, each = 5)
  rownames(v$fert_rat) <- 15:49
  v$cd4fert_rat <- projp$cd4fert_rat
  v$frr_art6mos <- projp$frr_art6mos
  v$frr_scalar <- projp$frr_scalar
  

  v
}
