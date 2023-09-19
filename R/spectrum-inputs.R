
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
  x <- c(6.663660,	6.627490,	6.587540,	6.550000,	6.515790,	6.480510,	6.442340,	6.399450,	6.350000,	6.292810,	6.226410,	6.148620,	6.057220,	5.950000,	5.795050,	5.582180,	5.341780,	5.104260,	4.900000,	4.728830,	4.568750,	4.417400,	4.272430,	4.131500,	3.988570,	3.844730,	3.708350,	3.587750,	3.491300,	3.414530,	3.346940,	3.287800,	3.236360,	3.191900,	3.150350,	3.110410,	3.075540,	3.049190,	3.034800,	3.029150,	3.025720,	3.023080,	3.019750,	3.014300,	3.002190,	2.981570,	2.955020,	2.925100,	2.894400,	2.859560,	2.817970,	2.773350,	2.729390,	2.689800,	2.654560,	2.621020,	2.588890,	2.557880,	2.527700,	2.498230,	2.469570)
  v$tfr = x


  adult_cd4_dist <- read.csv('tests/testdata/spectrum/v6.13/adult_cd4_dist.csv')
  adult_cd4_dist[is.na(adult_cd4_dist)] <- 0
  adult_cd4_dist <- adult_cd4_dist[,2:7]
  adult_cd4_dist_array <- array(unlist(adult_cd4_dist), dim = c(7,6), dimnames = list(cd4_adult = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'),
                                                                              cd4_adol = c('gte1000', '750-1000', '500-749', '350-499', '200-349', 'lte200')))

  v$adult_cd4_dist <- adult_cd4_dist_array


  v$mtct_trans <- (c(0.15, 0.15, 0.27, 0.27, 0.37, 0.37, 0.37))

  v$fert_mult_by_age <- rep(c(1.153260, 1.001870, 0.909590,0.912760, 0.883990, 0.883990, 0.883990), each = 5)

  v$fert_mult_offart <- c(1, 0.96, 0.88, 0.78, 0.61, 0.38, 0.30)

  v$fert_mult_onart <- rep(c(1.094460, 1.006870, 0.916170, 0.810910, 0.637220, 0.637220, 0.637220), each = 5)


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
  v$cd4fert_rat <- projp$cd4fert_rat
  v$frr_art6mos <- projp$frr_art6mos
  v$frr_scalar <- projp$frr_scalar
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



  v
}



prepare_hc_leapfrog_projp <- function(pjnz, params){
  projp <- read_hivproj_param(pjnz)
  ## Hard coded to expand age groups 15-24, 25-34, 35-44, 45+ to
  ## single-year ages 15:80.
  ## Requires extension for coarse HIV age group stratification
  idx_expand_full <- rep(1:4, times = c(10, 10, 10, 36))
  idx_expand_coarse <- rep(1:4, times = c(3, 2, 2, 2))

  v = params
  ## paed input
  v$paed_incid_input <- projp$nosocom_infections_04
  ## Hardcoded, this is putting all individuals in the highest cd4 category bc i think thats how the nosocomial infections work
  v$paed_cd4_dist <- c(0.6, 0.12, 0.1, 0.09, 0.05, 0.03, 0.01)
  ## v$paed_cd4_dist <- c(1, 0, 0, 0, 0, 0, 0)
  paed_cd4_prog <- array(data = 0, dim = c(7), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', 'lte5')) )
  paed_cd4_prog[] <- c(0.14, 0.37, 0.3, 0.35, 0.4, 0.4, 0)

  adol_cd4_prog <- array(data = 0, dim = c(6), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', 'lte200')) )
  adol_cd4_prog[] <- c(0.3028, 0.3028, 0.2575, 0.2122, 0.1669, 0)

  v$paed_cd4_prog <- paed_cd4_prog
  v$adol_cd4_prog <- adol_cd4_prog

  paed_cd4_mort <- array(data = 0, dim = c(7, 4, 5), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                     transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                     age = c(0:4)))
  ## 0-2
  paed_cd4_mort[,1,1:3] <- c(0.25801, 0.31513, 0.38490, 0.47011, 0.57418, 0.70129, 0.85655)
  paed_cd4_mort[,2,1:3] <- c(0.15385, 0.18791, 0.22951, 0.28031, 0.34237, 0.41817, 0.51074)
  paed_cd4_mort[,3,1:3] <- rep(0.08743, times = 7)
  paed_cd4_mort[,4,1:3] <- rep(0.02450, times = 7)

  ## 3-4
  paed_cd4_mort[,1,4:5] <- c(0.06470, 0.07902, 0.09652, 0.11789, 0.14398, 0.175860, 0.21479)
  paed_cd4_mort[,2,4:5] <- rep(0.04622, times = 7)
  paed_cd4_mort[,3,4:5] <- c(0.03383, 0.03485, 0.03609, 0.03761, 0.03946, 0.04172, 0.04449)
  paed_cd4_mort[,4,4:5] <- c(0.02144, 0.02347, 0.02596, 0.02899, 0.03270, 0.03723, 0.04276)





  paed_art_mort <- array(data = 0, dim = c(7, 3, 5), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                     transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                     age = c(0:4)))
  ## 0-6 mo on treatment
  paed_art_mort[,1,1] <- c(0.133460, 0.142490, 0.154150, 0.170110, 0.194100, 0.237140, 0.364790)
  paed_art_mort[,1,2:3] <- c(0.061010, 0.065140, 0.070470, 0.077770, 0.088730, 0.108400, 0.166760)
  paed_art_mort[,1,4:5] <- c(0.028720, 0.030660, 0.033170, 0.036600, 0.041760, 0.051020, 0.078490)

  ## 7-12 mo on treatment
  paed_art_mort[,2,1] <- c(0.048660, 0.051950, 0.056200, 0.062020, 0.070770, 0.086460, 0.133000)
  paed_art_mort[,2,2:3] <- c(0.022240, 0.023750, 0.025690, 0.028350, 0.032350, 0.039520, 0.060800)
  paed_art_mort[,2,4:5] <- c(0.010470, 0.011180, 0.012090, 0.013340, 0.015230, 0.018600, 0.028620)

  ## 12+ mo on treatment
  paed_art_mort[,3,1] <- rep(0, 7)
  paed_art_mort[,3,2:3] <- c(0.018610, 0.019340, 0.020260, 0.021480, 0.023220, 0.026140, 0.033720)
  paed_art_mort[,3,4:5] <- c(0.008760, 0.009100, 0.009540, 0.010110, 0.010930, 0.012300, 0.015870)

  adol_cd4_mort <- array(data = 0, dim = c(6, 4, 10), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                      transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                      age =5:14))
  ## 5 - 14
  adol_cd4_mort[,1,] <- c(0.02617, 0.03197, 0.03905, 0.04769, 0.05825, 0.07114)
  adol_cd4_mort[,2,] <- c(0.02572, 0.03142, 0.03837, 0.04687, 0.05724, 0.06991)
  adol_cd4_mort[,3,] <- c(0.02409, 0.02942, 0.03593, 0.04388, 0.05360, 0.06547)
  adol_cd4_mort[,4,] <- c(0.02245, 0.02742, 0.03349, 0.04090, 0.04996, 0.06102)

  adol_art_mort <- array(data = 0, dim = c(6, 3, 10), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                      transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                      age = c(5:14)))
  ## 0-6 mo on treatment
  adol_art_mort[,1,1:5] <- c(0.013800, 0.015860, 0.019110, 0.023650, 0.030090, 0.052640)
  adol_art_mort[,1,6:10] <- c(0.015140, 0.017400, 0.020950, 0.025930, 0.032990, 0.057730)

  ## 7-12 mo on treatment
  adol_art_mort[,2,1:5] <- c(0.005180, 0.005960, 0.007170, 0.008880, 0.011300, 0.019770)
  adol_art_mort[,2,6:10] <- c(0.005680, 0.006530, 0.007870, 0.009740, 0.012390, 0.021680)

  ## 12+ mo on treatment
  adol_art_mort[,3,1:5] <- c(0.005240, 0.005620, 0.006170, 0.006880, 0.007770, 0.010320)
  adol_art_mort[,3,6:10] <- c(0.005740, 0.006160, 0.006770, 0.007540, 0.008520, 0.011320)

  v$paed_cd4_mort <- paed_cd4_mort
  v$adol_cd4_mort <- adol_cd4_mort
  v$paed_art_mort <- paed_art_mort
  v$adol_art_mort <- adol_art_mort

  art_dist_paed <- read.csv('tests/testdata/spectrum/v6.13/paed_art_dist.csv')
  expand <- NULL
  for(yr in 2023:2030){
    x = art_dist_paed %>% dplyr::filter(year == 2022)
    x = x %>% dplyr::mutate(year = yr)
    expand <- rbind(expand, x)
  }
  art_dist_paed <- rbind(art_dist_paed, expand)
  art_dist_paed <- array(data = art_dist_paed$value, dim = c(15, length(unique(art_dist_paed$year))), dimnames = list( age= 0:14, year = sort(unique(art_dist_paed$year))))
  v$art_dist_paed <- art_dist_paed

  ## pull in cotrim coverage numbers
  ctx_pct <- T
  ## ctx_effect_notrt <- c(rep(0.33, 5), rep(0,5))
  ctx_effect <- 0.33
  v$ctx_val <- input_childart(pjnz)$ctx
  v$ctx_effect <- ctx_effect

  ## pull in ART coverage numbers
  v$artpaeds_isperc <- rep(T,  length(1970:2029))
  v$paed_art_val <- input_childart(pjnz)$child_art

  ##PMTCT
  pmtct_list <- input_pmtct(pjnz)
  pmtct_list <- pmtct_list[c(3,4,1,2,5:7),,]
  v$pmtct <- pmtct_list

  if(sum(pmtct_list[,,1]) == 0){
    v$pmtct_input_isperc = rep(F, length(1970:2030))
  }else{
    v$pmtct_input_isperc = rep(T, length(1970:2030))
  }

  ##PMTCT dropout
  v$pmtct_dropout <-input_pmtct_retained(pjnz)

  ##rates of MTCT
  noart <- read.csv('tests/testdata/spectrum/v6.13/mtct_notrt.csv')
  noart$cd4 <- factor(x = noart$cd4, levels = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'))
  noart <- noart %>% arrange(cd4)
  noart <- noart %>% spread(key = pmtct, value = perinatal)
  noart_per <- noart %>% filter(type == 'perinatal') %>% select(-type)
  rownames(noart_per) <- unique(noart_per$cd4)
  noart_per <- noart_per %>% select(-cd4)
  noart_bf <- noart %>% filter(type != 'perinatal') %>% select(-type)
  rownames(noart_bf) <- unique(noart_bf$cd4)
  noart_bf <- noart_bf %>% select(-cd4)
  noart <- list(noart_per, noart_bf)
  v$pmtct_mtct <- array(unlist(noart), dim = c(7,5,2), dimnames = list(cd4 = rownames(noart[[1]]), time = colnames(noart[[1]]), trans_type = c('perinatal', 'bf')))


  art <- read.csv('tests/testdata/spectrum/v6.13/mtct_trt.csv')
  art$cd4 <- factor(x = art$cd4, levels = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'))
  art <- art %>% arrange(cd4)
  art <- art %>% spread(key = pmtct, value = perinatal)
  art <- art %>% select(cd4, type, `ART before pregnancy`, `ART >4 weeks before delivery`, `ART <4 weeks before delivery`)
  art_per <- art %>% filter(type == 'perinatal') %>% select(-type)
  rownames(art_per) <- unique(art_per$cd4)
  art_per <- art_per %>% select(-cd4)
  art_bf <- art %>% filter(type != 'perinatal') %>% select(-type)
  rownames(art_bf) <- unique(art_bf$cd4)
  art_bf <- art_bf %>% select(-cd4)
  art <- list(art_per, art_bf)
  v$art_mtct <- array(unlist(art), dim = c(7,3,2), dimnames = list(cd4 = rownames(art[[1]]), time = colnames(art[[1]]), trans_type = c('perinatal', 'bf')))


  mort_rr_art <- read.csv('tests/testdata/spectrum/v6.13/mort_RR.csv', header = F)
  mort_rr_art_target <- array(NA, dim = c(3, 15, 61), dimnames = list(transmission = c('0to6mo', '7to12mo', '12+mo'), age = 0:14, year = 1970:2030))
  mort_rr_art_target[1:2, 1:5,] <- rep(unlist(mort_rr_art[1,]), each = 10)
  mort_rr_art_target[3, 1:5,] <- rep(unlist(mort_rr_art[2,]), each = 5)
  mort_rr_art_target[1:2, 6:15,] <- rep(unlist(mort_rr_art[3,]), each = 20)
  mort_rr_art_target[3, 6:15,] <- rep(unlist(mort_rr_art[4,]), each = 10)
  v$mort_art_rr <- mort_rr_art_target

  art_dist_paed <- read.csv('tests/testdata/spectrum/v6.13/paed_art_dist.csv')
  expand <- NULL
  for(yr in 2023:2030){
    x = art_dist_paed %>% dplyr::filter(year == 2022)
    x = x %>% dplyr::mutate(year = yr)
    expand <- rbind(expand, x)
  }
  art_dist_paed <- rbind(art_dist_paed, expand)
  art_dist_paed <- array(data = art_dist_paed$value, dim = c(15, length(unique(art_dist_paed$year))), dimnames = list( age= 0:14, year = sort(unique(art_dist_paed$year))))
  v$init_art_dist <- art_dist_paed

  ##BF duration
  v$bf_duration <- input_breastfeeding_dur(pjnz)

  ## ART eligibility age, doing this in years rather than months
  v$paed_art_elig_age <- c(rep(0, 37), rep(1, 3), rep(2, 20))
  paed_art_elig_cd4 <- array(data = NA, dim = c(length(0:14), length(1970:2029)), dimnames = list(age = c(0:14), year = c(1970:2029)))
  ## corresponds with less than 25
  paed_art_elig_cd4[1:5,] <- 3
  ## correspond with 200 then 300
  paed_art_elig_cd4[6:15,1:40] <- 6
  paed_art_elig_cd4[6:15,41:60] <- 5
  v$paed_art_elig_cd4 <- paed_art_elig_cd4

  v$paed_art_ltfu <- rep(0, length(1970:2030))

  paed_cd4_transition <- array(0, dim = c(6,7), dimnames = list(cd4_count = c('gte1000', '750-1000', '500-749', '350-499', '200-349', 'lte200'), cd4_pct = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5')))
  paed_cd4_transition[1:6,1] <- c(0.608439, 0.185181, 0.105789, 0.055594, 0.018498, 0.026497)
  paed_cd4_transition[1:6,2] <- c(0.338733873, 0.222622262, 0.293529353, 0.093509351, 0.03550355, 0.01610161)
  paed_cd4_transition[1:6,3] <- c(0.2004, 0.2562, 0.3636, 0.1074, 0.0579, 0.0145)
  paed_cd4_transition[1:6,4] <- c(0.095, 0.1693, 0.3082, 0.2497, 0.1449, 0.0329)
  paed_cd4_transition[1:6,5] <- c(0.03880388, 0.090309031, 0.275927593, 0.259925993, 0.255725573, 0.079307931)
  paed_cd4_transition[1:6,6] <- c(0.018615705, 0.018615705, 0.099217744, 0.165065848, 0.363501337, 0.334983662)
  paed_cd4_transition[1:6,7] <- c(0, 0.0014, 0.00990099, 0.00710071, 0.04960496, 0.931993199)
  v$paed_cd4_transition <- paed_cd4_transition

  return(v)
}
