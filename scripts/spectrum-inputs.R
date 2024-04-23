
read_sx <- function(pjnz, use_ep5=FALSE) {

  if(use_ep5) {
    dpfile <- grep(".ep5$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  } else {
    dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  }

  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)

  exists_dptag <- function(tag, tagcol=1) {
    tag %in% dp[,tagcol]
  }

  dpsub <- function(tag, rows, cols, tagcol=1) {
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

read_netmigr <- function(pjnz, use_ep5=FALSE, adjust_u5mig = TRUE, sx = NULL) {

  if(use_ep5) {
    dpfile <- grep(".ep5$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  } else {
    dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  }

  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)

  exists_dptag <- function(tag, tagcol=1) {
    tag %in% dp[,tagcol]
  }
  dpsub <- function(tag, rows, cols, tagcol=1) {
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

  demp$births_sex_prop <- rbind(male = demp$srb, female = 100) / (demp$srb + 100)

  ## normalise ASFR distribution
  demp$asfr <- sweep(demp$asfr, 2, demp$tfr / colSums(demp$asfr), "*")


  ## NOTE: Reading this to obtain the Spectrum version number
  ##       This is a lot of redundant effort.
  projp <- eppasm::read_hivproj_param(pjnz)
  if (!grepl("^[4-6]\\.[0-9]", projp$spectrum_version)) {
    stop(paste0("Spectrum version not recognized: ", projp$spectrum_version))
  }
  demp$projection_period <- if (projp$spectrum_version >= "6.2") {
    "calendar"
  } else {
    "midyear"
  }

  if (demp$projection_period == "midyear") {
    demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  } else {
    demp$netmigr_adj <- demp$netmigr
  }

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

  projp <- eppasm::read_hivproj_param(pjnz)

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
  # x <- c(6.663660,	6.627490,	6.587540,	6.550000,	6.515790,	6.480510,	6.442340,	6.399450,	6.350000,	6.292810,	6.226410,	6.148620,	6.057220,	5.950000,	5.795050,	5.582180,	5.341780,	5.104260,	4.900000,	4.728830,	4.568750,	4.417400,	4.272430,	4.131500,	3.988570,	3.844730,	3.708350,	3.587750,	3.491300,	3.414530,	3.346940,	3.287800,	3.236360,	3.191900,	3.150350,	3.110410,	3.075540,	3.049190,	3.034800,	3.029150,	3.025720,	3.023080,	3.019750,	3.014300,	3.002190,	2.981570,	2.955020,	2.925100,	2.894400,	2.859560,	2.817970,	2.773350,	2.729390,	2.689800,	2.654560,	2.621020,	2.588890,	2.557880,	2.527700,	2.498230,	2.469570)
  # v$tfr = x


  adult_cd4_dist <- array(data = 0, dim = c(7,6), dimnames = list(adult_cd4_categories = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99','<50'),
                                                                   hc2_cd4_categories = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200')))
  adult_cd4_dist[1,1:3] <- 1
  adult_cd4_dist[2,4] <- 1
  adult_cd4_dist[3:4,5] <- c(0.6665589, 1-0.6665589)
  adult_cd4_dist[5:7,6] <- c(0.35, 0.21, 0.44)


  adult_cd4_dist_array <-adult_cd4_dist

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


#' @export
prepare_hc_leapfrog_projp <- function(pjnz, params, pop_1){
  ## Hard coded to expand age groups 15-24, 25-34, 35-44, 45+ to
  ## single-year ages 15:80.
  ## Requires extension for coarse HIV age group stratification
  idx_expand_full <- rep(1:4, times = c(10, 10, 10, 36))
  idx_expand_coarse <- rep(1:4, times = c(3, 2, 2, 2))

  v = params
  ## paed input
  v$paed_incid_input <- leapfrog:::dp_read_nosocom_infections(pjnz)
  v$paed_cd4_dist <- leapfrog:::dp_read_paed_cd4_dist(pjnz) / 100

  prog = leapfrog:::dp_read_paed_cd4_prog(pjnz)
  paed_cd4_prog <- array(c(prog[1,1:6],0), dim = 7, dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5')))

  adol_cd4_prog <- array(c(prog[1,14:18],0), dim = 6, dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', 'lte200')))

  v$paed_cd4_prog <- paed_cd4_prog
  v$adol_cd4_prog <- adol_cd4_prog

  mort <- leapfrog:::dp_read_paed_cd4_mort(pjnz)
  paed_cd4_mort <- array(data = 0, dim = c(7, 4, 5), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                     transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                     age = c(0:4)))
  ## 0-2
  paed_cd4_mort[,1,1:3] <- mort[1,]
  paed_cd4_mort[,2,1:3] <- mort[2,]
  paed_cd4_mort[,3,1:3] <- mort[3,]
  paed_cd4_mort[,4,1:3] <- mort[4,]

  ## 3-4
  paed_cd4_mort[,1,4:5] <- mort[5,]
  paed_cd4_mort[,2,4:5] <- mort[6,]
  paed_cd4_mort[,3,4:5] <- mort[7,]
  paed_cd4_mort[,4,4:5] <- mort[8,]

  adol_cd4_mort <- array(data = 0, dim = c(6, 4, 10), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                      transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                      age =5:14))
  ## 5 - 14
  adol_cd4_mort[,1,] <- mort[9,2:7]
  adol_cd4_mort[,2,] <- mort[10,2:7]
  adol_cd4_mort[,3,] <- mort[11,2:7]
  adol_cd4_mort[,4,] <- mort[12,2:7]


  mort <- leapfrog:::dp_read_paed_art_mort(pjnz)
  paed_art_mort <- array(data = 0, dim = c(7, 3, 5), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                     transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                     age = c(0:4)))
  ## 0-6 mo on treatment
  paed_art_mort[,1,1] <- mort$hc1_lt6[,1,1]
  paed_art_mort[,1,2:3] <- mort$hc1_lt6[,2,1]
  paed_art_mort[,1,4:5] <- mort$hc1_lt6[,3,1]

  ## 7-12 mo on treatment
  paed_art_mort[,2,1] <- mort$hc1_6to12[,1,1]
  paed_art_mort[,2,2:3] <- mort$hc1_6to12[,2,1]
  paed_art_mort[,2,4:5] <- mort$hc1_6to12[,3,1]

  ## 12+ mo on treatment
  paed_art_mort[,3,1] <- mort$hc1_gte12[,1,1]
  paed_art_mort[,3,2:3] <- mort$hc1_gte12[,2,1]
  paed_art_mort[,3,4:5] <- mort$hc1_gte12[,3,1]


  adol_art_mort <- array(data = 0, dim = c(6, 3, 10), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                      transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                      age = c(5:14)))
  ## 0-6 mo on treatment
  adol_art_mort[,1,1:5] <- mort$hc2_lt6[,1,1]
  adol_art_mort[,1,6:10] <- mort$hc2_lt6[,2,1]

  ## 7-12 mo on treatment
  adol_art_mort[,2,1:5] <- mort$hc2_6to12[,1,1]
  adol_art_mort[,2,6:10] <- mort$hc2_6to12[,2,1]

  ## 12+ mo on treatment
  adol_art_mort[,3,1:5] <- mort$hc2_gte12[,1,1]
  adol_art_mort[,3,6:10] <- mort$hc2_gte12[,2,1]

  v$paed_cd4_mort <- paed_cd4_mort
  v$adol_cd4_mort <- adol_cd4_mort
  v$paed_art_mort <- paed_art_mort
  v$adol_art_mort <- adol_art_mort

  mtct_rates_input <- leapfrog:::dp_read_mtct_rates(pjnz)
  art_mtct <- array(0, dim = c(7,3,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'),
                                                       time = c('ART <4 weeks before delivery', 'ART >4 weeks before delivery', 'ART before pregnancy'), trans_type = c('perinatal', 'bf')))
  art_mtct[,1,1] <- mtct_rates_input[11,1] / 100
  art_mtct[,2,1] <- mtct_rates_input[10,1] / 100
  art_mtct[,3,1] <- mtct_rates_input[9,1] / 100
  art_mtct[5:7,1,2] <- mtct_rates_input[11,2] / 100
  art_mtct[5:7,2,2] <- mtct_rates_input[10,2] / 100
  art_mtct[5:7,3,2] <- mtct_rates_input[9,2] / 100
  v$art_mtct <- art_mtct

  art_dist_paed <- leapfrog:::dp_read_art_dist(pjnz)
  v$art_dist_paed <- art_dist_paed

  ## pull in cotrim coverage numbers
  ctx_pct <- leapfrog:::input_childart(pjnz)$ctx_percent
  ctx_pct[is.na(ctx_pct)] <- FALSE
  v$ctx_val_ispercent <- ctx_pct
  ## ctx_effect_notrt <- c(rep(0.33, 5), rep(0,5))
  ctx_effect <- rep(0.33,61)
  v$ctx_val <- leapfrog:::input_childart(pjnz)$ctx
  if(any(v$ctx_val_ispercent)){
    v$ctx_val[v$ctx_val_ispercent] <- v$ctx_val[v$ctx_val_ispercent] / 100
  }
  v$ctx_effect <- as.array(ctx_effect)

  ## pull in ART coverage numbers
  art = input_childart(pjnz)
  v$artpaeds_isperc <- art$art_ispercent[2,]
  v$artpaeds_isperc[] <- as.integer(ifelse(art$art_ispercent[2,] == FALSE, 0, 1))
  art$child_art[1,which(art$art_ispercent[2,])] <- art$child_art[1,which(art$art_ispercent[2,])] / 100
  v$paed_art_val <- art$child_art
  v$paed_art_age_spec <- art$age_spec

  ##PMTCT
  pmtct_list <- leapfrog:::input_pmtct(pjnz)
  pmtct_list <- pmtct_list[c(3,4,1,2,5:7),,]
  v$pmtct <- pmtct_list

  if(sum(pmtct_list[,,1]) == 0){
    v$pmtct_input_isperc = rep(F, length(1970:2030))
  }else{
    v$pmtct_input_isperc = rep(T, length(1970:2030))
  }

  v$pmtct_input_isperc <- !(apply(leapfrog:::input_pmtct_ispercent(pjnz), 2, any))

  ##PMTCT dropout
  v$pmtct_dropout <- leapfrog:::input_pmtct_retained(pjnz)

  ##rates of MTCT
  mtct_trt <- array(data = 0, dim = c(7,7,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'),
                                                       pmtct_reg = c('option A', 'option B', 'single dose nevirapine', 'WHO 2006 dual ARV regimen', 'ART before pregnancy',
                                                                     'ART >4 weeks before delivery', 'ART <4 weeks before delivery'),
                                                       transmission_type = c('perinatal', 'breastfeeding')))
  mtct_trt[,1,1] <-  mtct_rates_input[7,1] /100
  mtct_trt[,2,1] <-  mtct_rates_input[8,1]/100
  mtct_trt[,3,1] <- mtct_rates_input[5,1]/100
  mtct_trt[,4,1] <- mtct_rates_input[6,1]/100
  mtct_trt[5:7,1,2] <- mtct_rates_input[7,3]/100
  mtct_trt[5:7,2,2] <- mtct_rates_input[8,3]/100
  mtct_trt[1:4,3,2] <- mtct_rates_input[5,2]/100
  mtct_trt[5:7,3,2] <- mtct_rates_input[5,3]/100
  mtct_trt[,4,2] <- mtct_rates_input[6,3]/100
  mtct_trt[,5,1] <- mtct_rates_input[9,1]/100
  mtct_trt[,6,1] <- mtct_rates_input[10,1]/100
  mtct_trt[,7,1] <- mtct_rates_input[11,1]/100
  mtct_trt[5:7,5,2] <- mtct_rates_input[9,2]/100
  mtct_trt[5:7,6,2] <- mtct_rates_input[10,2]/100
  mtct_trt[5:7,7,2] <- mtct_rates_input[11,2]/100
  v$pmtct_mtct <- mtct_trt

  mtct <- array(data = NA, dim = c(7,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'), trans_type = c('perinatal', 'bf')))
  mtct[,1] <- c(mtct_rates_input[3,1], mtct_rates_input[3,1],
                mtct_rates_input[2,1], mtct_rates_input[2,1],
                rep(mtct_rates_input[1,1],3)) /100
  mtct[,2] <- c(rep(mtct_rates_input[3,3],2), rep(mtct_rates_input[2,2],2), rep(mtct_rates_input[1,2],3)) / 100
  v$mtct <- mtct

  mort_rr_art <- leapfrog:::dp_read_child_mort_mult(pjnz)
  mort_rr_art_target <- array(NA, dim = c(3, 15, 61), dimnames = list(transmission = c('0to6mo', '7to12mo', '12+mo'), age = 0:14, year = 1970:2030))
  mort_rr_art_target[1:2, 1:5,] <- rep(unlist(mort_rr_art[1,]), each = 10)
  mort_rr_art_target[3, 1:5,] <- rep(unlist(mort_rr_art[2,]), each = 5)
  mort_rr_art_target[1:2, 6:15,] <- rep(unlist(mort_rr_art[3,]), each = 20)
  mort_rr_art_target[3, 6:15,] <- rep(unlist(mort_rr_art[4,]), each = 10)
  v$mort_art_rr <- mort_rr_art_target

  art_dist_paed <- leapfrog:::dp_read_art_dist(pjnz)
  v$init_art_dist <- art_dist_paed

  ##BF duration
  bf_duration <- leapfrog:::input_breastfeeding_dur(pjnz)
  v$bf_duration_art <- bf_duration[,,2]
  v$bf_duration_no_art <- bf_duration[,,1]
  ##only keeping this for leapfrog
  v$bf_duration = bf_duration

  art_elig = leapfrog:::dp_read_paed_art_eligibility(pjnz)
  v$paed_art_elig_age <- art_elig$age_elig / 12 ##converts from months to years

  ##MKW: stopped here
  cd4_elig <- art_elig$cd4_elig[c(5:7,4),]
  ##Changing the input from CD4 count or percentages to ordinal categories
  ###Easier to do it here than in the leapfrog code
  get_ordinal <- function(vec){
    out <- c()
    for (i in 1:length(vec)) {
      if(vec[i] >= 1000){
        val = 1
      }else if( vec[i] >= 750){
        val = 2
      }else if(vec[i] >= 500 ){
        val = 3
      }else if(vec[i] >= 350 ){
        val = 4
      }else if(vec[i] >= 200 ){
        val = 5
      }else if(vec[i] >= 50 ){
        val = 6
      }else if(vec[i] > 30 ){
        val = 0
      }else if(vec[i] >= 26 ){
        val = 1
      }else if(vec[i] >= 21 ){
        val = 2
      }else if(vec[i] >= 16 ){
        val = 3
      }else if(vec[i] >= 11 ){
        val = 4
      }else if(vec[i] >= 5 ){
        val = 5
      }else {
        val = 6
      }
      out[i] <- val
    }
    out
  }

  paed_art_elig_cd4 <- array(data = NA, dim = c(15, 61), dimnames = list(age = c(0:14), year = c(1970:2030)))
  paed_art_elig_cd4[1,] <- get_ordinal(unname(cd4_elig[1,]))
  paed_art_elig_cd4[2:3,] <- get_ordinal(unname(cd4_elig[2,]))
  paed_art_elig_cd4[4:5,] <- get_ordinal(unname(cd4_elig[3,]))
  paed_art_elig_cd4[6:15,] <- rep(get_ordinal(cd4_elig[4,]), each = length(6:15))
  v$paed_art_elig_cd4 <- paed_art_elig_cd4


  v$paed_art_ltfu <- leapfrog:::input_childart_ltfu(pjnz) / 100

  paed_cd4_transition <- array(0, dim = c(6,7), dimnames = list(cd4_count = c('gte1000', '750-1000', '500-749', '350-499', '200-349', 'lte200'), cd4_pct = c('gte30', '26-30', '21-25', '16-20', '11-15', '5-10', 'lte5')))
  paed_cd4_transition[1:6,1] <- c(0.608439, 0.185181, 0.105789, 0.055594, 0.018498, 0.026497)
  paed_cd4_transition[1:6,2] <- c(0.338733873, 0.222622262, 0.293529353, 0.093509351, 0.03550355, 0.01610161)
  paed_cd4_transition[1:6,3] <- c(0.2004, 0.2562, 0.3636, 0.1074, 0.0579, 0.0145)
  paed_cd4_transition[1:6,4] <- c(0.095, 0.1693, 0.3082, 0.2497, 0.1449, 0.0329)
  paed_cd4_transition[1:6,5] <- c(0.03880388, 0.090309031, 0.275927593, 0.259925993, 0.255725573, 0.079307931)
  paed_cd4_transition[1:6,6] <- c(0.018615705, 0.018615705, 0.099217744, 0.165065848, 0.363501337, 0.334983662)
  paed_cd4_transition[1:6,7] <- c(0, 0.0014, 0.00990099, 0.00710071, 0.04960496, 0.931993199)
  paed_cd4_transition[6,1] <- 1 - sum(paed_cd4_transition[-6,1])
  paed_cd4_transition[6,7] <- 1 - sum(paed_cd4_transition[-6,7])
  # paed_cd4_transition[1:2,1] <- c(0.71, 0.29)
  # paed_cd4_transition[2:3,2] <- c(0.6, 0.4)
  # paed_cd4_transition[3:4,3] <- c(0.83, 0.17)
  # paed_cd4_transition[4:5,4] <- c(0.77, 0.23)
  # paed_cd4_transition[5:6,5] <- c(0.89, 0.11)
  # paed_cd4_transition[6,6] <- c(1)
  # paed_cd4_transition[6,7] <- c(1)

  v$paed_cd4_transition <- paed_cd4_transition


  ## projection parameters
  dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }
  yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1

  abort <- dpsub(tag = "<PregTermAbortionPerNum MV2>", rows = 2, cols = timedat.idx)
  abort <- as.numeric(abort) / 100 / 100
  v$abortions <- abort


  ##extract needed outputs to just run paed model
  wlhiv_births <- dpsub("<ChildNeedPMTCT MV>", 2, timedat.idx) %>% unname()
  names(wlhiv_births) <- proj.years
  rownames(wlhiv_births) <- NULL
  v$mat_hiv_births <- as.array(as.numeric(unlist(wlhiv_births)))

  pjnz1 <- pjnz
  specres <- eppasm::read_hivproj_output(pjnz1)
  newinf <- specres$newinf.f[4:10,] %>% colSums()
  newinf_rate <- newinf / colSums(specres$totpop.f[1:10,])
  v$incrate <- as.array(as.numeric(unlist(newinf_rate)))

  fp1 <- eppasm::prepare_directincid(pjnz1)
  fp1$tARTstart <- 61L
  #Need to run the adult model to pull out the propotion by CD4 category
  mod1 <- eppasm::simmod(fp1)

  #cd4
  wlhiv_cd4 <- array(as.numeric(unlist(dpsub("<CD4Distribution15_49 MV2>", 19:25, timedat.idx))), dim = c(7,length(timedat.idx)))
  v$prop_gte350 <- colSums(wlhiv_cd4[1:2,]) / colSums(wlhiv_cd4)
  v$prop_lt200 <- colSums(wlhiv_cd4[5:7,]) / colSums(wlhiv_cd4)

  v$hc_age_coarse <- rep(c(1,2,3), each = 5)


  return(v)
}
