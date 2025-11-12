#' Read Spectrum .DP file
#'
#' This function returns the Spectrum .DP file read as character CSV. Not intended
#' for direct use, but passing to other functions to parse.
#'
#' @param pjnz file path to Spectrum PJNZ file
#'
#' @return Matrix with class "spectrum_dp". Not intended for direct use.
#'
#' @examples
#'
#' pjnz <- system.file("pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ", package = "frogger")
#' dp <- read_dp(pjnz)
#' class(dp)
#' @noRd
read_dp <- function(pjnz) {

  stopifnot(grepl("\\.(pjnz|zip)$", pjnz, ignore.case = TRUE))

  dpfile <- grep("\\.DP$", utils::unzip(pjnz, list = TRUE)$Name, value = TRUE)
  stopifnot(length(dpfile) == 1)

  dp <- utils::read.csv(unz(pjnz, dpfile), as.is = TRUE)
  class(dp) <- c(class(dp), "spectrum_dp")

  dp
}


## TODO: figure out how to do roxygen documentation
prepare_abortion_input <- function(data){
  yr_start <- data$first_year$data
  yr_end <- data$final_year$data
  proj.years <- yr_start:yr_end

  abortion <- array(0, dim = c(2, length(proj.years)), dimnames = list(value = c('Input', 'Percent'), year = proj.years))
  abortion["Input",] <- data$preg_term_abortion$data
  abortion["Percent",] <- data$preg_term_abortion_pernum$data

  return(abortion)
}

prepare_cd4_progression <- function(data){
  prog = data$child_ann_rate_progress_lower_cd4$data
  hc1_cd4_prog <- array(0, dim = c(7,2,2), dimnames = list(cd4pct = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                            age = c('0-2', '3-4'),
                                                            sex = c('male', 'female')))
  hc2_cd4_prog <- array(0, dim = c(6,1,2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', 'lte200'),
                                                            age = c('5-14'),
                                                            sex = c('male', 'female')))
  hc1_cd4_prog[,"0-2","male"] <- c(prog["male",grepl("Age: 0-2,", colnames(prog))], 0)
  hc1_cd4_prog[,"0-2","female"] <- c(prog["female",grepl("Age: 0-2,", colnames(prog))], 0)
  hc1_cd4_prog[,"3-4","male"] <- c(prog["male",grepl("Age: 3-4,", colnames(prog))], 0)
  hc1_cd4_prog[,"3-4","female"] <- c(prog["female",grepl("Age: 3-4,", colnames(prog))], 0)
  hc2_cd4_prog[,"5-14","male"] <- c(prog["male",grepl("Age: 5-14,", colnames(prog))], 0)
  hc2_cd4_prog[,"5-14","female"] <- c(prog["female",grepl("Age: 5-14,", colnames(prog))], 0)

  return(list(hc1_cd4_prog = hc1_cd4_prog,
              hc2_cd4_prog = hc2_cd4_prog))
}

prepare_no_art_mort <- function(data){
  mort <- data$child_mort_by_cd4_no_art$data
  hc1_cd4_mort <- array(data = 0, dim = c(7, 4, 5), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                     transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                     age = c(0:4)))
  hc2_cd4_mort <- array(data = 0, dim = c(6, 4, 10), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                      transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+'),
                                                                      age =5:14))
  ## 0-2
  hc1_cd4_mort[,"perinatal",c("0", "1", "2")] <- mort[grepl("Age 0-2: Perinatal", rownames(mort)),]
  hc1_cd4_mort[,"bf0-6",c("0", "1", "2")] <- mort[grepl("Age 0-2: Breastfeeding, <6MOS", rownames(mort)),]
  hc1_cd4_mort[,"bf7-12",c("0", "1", "2")] <- mort[grepl("Age 0-2: Breastfeeding, 7-12MOS", rownames(mort)),]
  hc1_cd4_mort[,"bf12+",c("0", "1", "2")] <- mort[grepl("Age 0-2: Breastfeeding, >12MOS", rownames(mort)),]

  ## 3-4
  hc1_cd4_mort[,"perinatal",c("3", "4")] <- mort[grepl("Age 3-4: Perinatal", rownames(mort)),]
  hc1_cd4_mort[,"bf0-6",c("3", "4")] <- mort[grepl("Age 3-4: Breastfeeding, <6MOS", rownames(mort)),]
  hc1_cd4_mort[,"bf7-12",c("3", "4")] <- mort[grepl("Age 3-4: Breastfeeding, 7-12MOS", rownames(mort)),]
  hc1_cd4_mort[,"bf12+",c("3", "4")] <- mort[grepl("Age 3-4: Breastfeeding, >12MOS", rownames(mort)),]

  ## 5-14, skip first index to account for the difference in CD4 categories
  hc2_cd4_mort[,"perinatal",as.character(5:14)] <- mort[grepl("Age 5-14: Perinatal", rownames(mort)),][-1]
  hc2_cd4_mort[,"bf0-6",as.character(5:14)] <- mort[grepl("Age 5-14: Breastfeeding, <6MOS", rownames(mort)),][-1]
  hc2_cd4_mort[,"bf7-12",as.character(5:14)] <- mort[grepl("Age 5-14: Breastfeeding, 7-12MOS", rownames(mort)),][-1]
  hc2_cd4_mort[,"bf12+",as.character(5:14)] <- mort[grepl("Age 5-14: Breastfeeding, >12MOS", rownames(mort)),][-1]

  return(list(hc1_cd4_mort = hc1_cd4_mort,
              hc2_cd4_mort = hc2_cd4_mort))
}

prepare_art_mort <- function(data){
  mort_lt6 <- data$child_mort_by_cd4_with_art_0to6$data
  mort_7to12 <- data$child_mort_by_cd4_with_art_7to12$data
  mort_gt12 <- data$child_mort_by_cd4_with_art_gt12$data
  ##TODO: stratify this by sex
  hc1_art_mort <- array(data = 0, dim = c(7, 3, 5), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                                     transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                     age = c(0:4)))
  hc2_art_mort <- array(data = 0, dim = c(6, 3, 10), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'),
                                                                      transmission = c('0to6mo', '7to12mo', '12+mo'),
                                                                      age = c(5:14)))
  ## 0-6 mo on treatment
  hc1_art_mort[,'0to6mo',"0"] <- mort_lt6["male",grepl("0: ", colnames(mort_lt6))]
  hc1_art_mort[,'0to6mo',as.character(1:2)] <- mort_lt6["male",grepl("1-2: ", colnames(mort_lt6))]
  hc1_art_mort[,'0to6mo',as.character(3:4)] <- mort_lt6["male",grepl("3-4: ", colnames(mort_lt6))]
  hc2_art_mort[,"0to6mo",as.character(5:9)] <- mort_lt6["male",grepl("5-9: ", colnames(mort_lt6))]
  hc2_art_mort[,"0to6mo",as.character(10:14)] <- mort_lt6["male",grepl("10-14: ", colnames(mort_lt6))]

  ## 7-12 mo on treatment
  hc1_art_mort[,"7to12mo","0"] <- mort_7to12["male",grepl("0: ", colnames(mort_7to12))]
  hc1_art_mort[,"7to12mo",as.character(1:2)] <- mort_7to12["male",grepl("1-2: ", colnames(mort_7to12))]
  hc1_art_mort[,"7to12mo",as.character(3:4)] <- mort_7to12["male",grepl("3-4: ", colnames(mort_7to12))]
  hc2_art_mort[,"7to12mo",as.character(5:9)] <- mort_7to12["male",grepl("5-9: ", colnames(mort_7to12))]
  hc2_art_mort[,"7to12mo",as.character(10:14)] <- mort_7to12["male",grepl("10-14: ", colnames(mort_7to12))]

  ## 12+ mo on treatment
  hc1_art_mort[,"12+mo","0"] <- mort_gt12["male",grepl("0: ", colnames(mort_gt12))]
  hc1_art_mort[,"12+mo",as.character(1:2)] <- mort_gt12["male",grepl("1-2: ", colnames(mort_gt12))]
  hc1_art_mort[,"12+mo",as.character(3:4)] <- mort_gt12["male",grepl("3-4: ", colnames(mort_gt12))]
  hc2_art_mort[,"12+mo",as.character(5:9)] <- mort_gt12["male",grepl("5-9: ", colnames(mort_gt12))]
  hc2_art_mort[,"12+mo",as.character(10:14)] <- mort_gt12["male",grepl("10-14: ", colnames(mort_gt12))]

  return(list(hc1_art_mort = hc1_art_mort,
              hc2_art_mort = hc2_art_mort))
}

prepare_cotrim_effect <- function(data){
  ##cotrim is effective for five years for children not on ART and for four years for children on ART
  ctx_effect <- data$effect_treat_child$data
  off_art_ctx <- sum(as.numeric(unlist(ctx_effect["no art",]))) / 5
  on_art_ctx.lte12mo <- sum(as.numeric(unlist(ctx_effect["art",1])))
  on_art_ctx.gte12mo <- sum(as.numeric(unlist(ctx_effect["art",2:5]))) / 4
  ctx_effect <- array(data = c(off_art_ctx, on_art_ctx.lte12mo, on_art_ctx.gte12mo),
                      dim = c(3),
                      dimnames = list(ctx_effect = c('Off ART', 'On ART, lte12mo', 'On ART, gte12mo')))
  return(ctx_effect)
}

prepare_pmtct <- function(data){
  pmtct <- data$arv_regimen$data
  pmtct[is.na(pmtct)] = 0
  pmtct_number <- pmtct[grepl('- Number', rownames(pmtct)) &
                          !grepl('postnatal', rownames(pmtct)) & !grepl('Postnatal', rownames(pmtct)) &
                          !grepl('Total', rownames(pmtct)),]
  pmtct_pct <- pmtct[grepl('- Percent', rownames(pmtct)) &
                       !grepl('postnatal', rownames(pmtct)) & !grepl('Postnatal', rownames(pmtct)) &
                       !grepl('Total', rownames(pmtct)),]
  pmtct_pct <- pmtct_pct[-(which(rownames(pmtct_pct) == "No prophylaxis- Percent")),]
  pmtct_input_isperc <- rep(1, ncol(pmtct_number))
  pmtct_input_isperc[colSums(pmtct_number) > 0] <- 0
  order =  c("Option A", "Option B", "Single dose nevirapine",
             "WHO 2006 dual ARV regimen", "ART: Started before pregnancy",
             "ART: Started during pregnancy >4 weeks", "ART: Started during pregnancy <4 weeks")
  rownames(pmtct_number) <- gsub(pattern = "- Number", replacement = "", rownames(pmtct_number))
  rownames(pmtct_pct) <- gsub(pattern = "- Percent", replacement = "", rownames(pmtct_pct))
  ##Order in the expected order for leapfrog
  pmtct_number <- pmtct_number[match(order, rownames(pmtct_number)),]
  pmtct_pct <- pmtct_pct[match(order, rownames(pmtct_pct)),]

  pmtct_new <- array(0, dim = c(7, 61), dimnames = list(pmtct = c("Option A", "Option B", "SDNVP", "Dual ARV", "Option B+: before pregnancy", "Option B+: >4 weeks", "Option B+: <4 weeks")))
  ## pick out which ones were inserted as numbers
  pmtct_new <- pmtct_number
  ## pick out which ones were inserted as percent
  pmtct_new[,which(pmtct_input_isperc == 1)] <- pmtct_pct[,which(pmtct_input_isperc == 1)]

  return(list(pmtct_new = pmtct_new,
              pmtct_input_isperc = pmtct_input_isperc))
}

prepare_pmtct_dropout <- function(data){
  yr_start <- data$first_year$data
  yr_end <- data$final_year$data
  proj.years <- yr_start:yr_end

  pmtct_dropout <- array(0, dim = c(7, 3, length(proj.years)),
                         dimnames = list(pmtct = c("Option A", "Option B", "SDNVP", "Dual ARV", "Option B+: before pregnancy", "Option B+: >4 weeks", "Option B+: <4 weeks"),
                                         drop_out_by = c('Delivery', '<12MOS breastfeeding', '>12MOS breastfeeding'),
                                         year = proj.years))
  pmtct_dropout[,"Delivery",] <- 100
  pmtct_dropout["Option B+: before pregnancy","Delivery",] <- data$percent_art_delivery$data["Percent already on ART retained at delivery",]
  pmtct_dropout["Option B+: >4 weeks","Delivery",] <- data$percent_art_delivery$data["Percent starting ART retained at delivery",]
  pmtct_dropout[c("Option A",
                  "Option B",
                  "Option B+: before pregnancy",
                  "Option B+: >4 weeks",
                  "Option B+: <4 weeks"),"<12MOS breastfeeding",] <- rep(data$arv_regimen$data[rownames(data$arv_regimen$data) %in%
                                                                                             c("Monthly dropout breastfeeding: ART 0-12 months breastfeeding" ),], each = 5)
  pmtct_dropout[c("Option A",
                  "Option B",
                  "Option B+: before pregnancy",
                  "Option B+: >4 weeks",
                  "Option B+: <4 weeks"),">12MOS breastfeeding",] <- rep(data$arv_regimen$data[rownames(data$arv_regimen$data) %in%
                                                                                             c("Monthly dropout breastfeeding: ART 12+ months breastfeeding" ),], each = 5)
  pmtct_dropout[is.na(pmtct_dropout)] <- 0
  pmtct_dropout <- pmtct_dropout / 100
  return(pmtct_dropout)
}

prepare_vertical_transmission <- function(data){
  order =  c("Option A", "Option B", "Single dose nevirapine",
             "WHO 2006 dual ARV regimen", "ART: Started before pregnancy",
             "ART: Started during pregnancy >4 weeks", "ART: Started during pregnancy <4 weeks")

  mtct <- data$trans_eff_assump$data / 100
  untrt_mtct <- mtct[!rownames(mtct) %in% order,]
  mtct <- mtct[match(order, rownames(mtct)),]
  mtct_trt <- array(data = 0, dim = c(7,7,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'),
                                                              pmtct_reg = c("Option A", "Option B", "SDNVP", "Dual ARV",
                                                                            "Option B+: before pregnancy", "Option B+: >4 weeks", "Option B+: <4 weeks"),
                                                              transmission_type = c('perinatal', 'breastfeeding')))
  mtct_trt[,"Option A","perinatal"] <- mtct["Option A", "Perinatal"]
  mtct_trt[,"Option B","perinatal"] <- mtct["Option B", "Perinatal"]
  mtct_trt[,"SDNVP","perinatal"] <- mtct["Single dose nevirapine", "Perinatal"]
  mtct_trt[,"Dual ARV","perinatal"] <- mtct["WHO 2006 dual ARV regimen", "Perinatal"]
  mtct_trt[,"Option B+: before pregnancy","perinatal"] <- mtct["ART: Started before pregnancy", "Perinatal"]
  mtct_trt[,"Option B+: >4 weeks","perinatal"] <- mtct["ART: Started during pregnancy >4 weeks", "Perinatal"]
  mtct_trt[,"Option B+: <4 weeks","perinatal"] <- mtct["ART: Started during pregnancy <4 weeks", "Perinatal"]

  mtct_trt[c('250-349', '200-249', '100-199', '50-99', '<50'),"Option A","breastfeeding"] <- mtct["Option A", "Breastfeeding (per month) >=350"]
  mtct_trt[c('250-349', '200-249', '100-199', '50-99', '<50'),"Option B","breastfeeding"] <- mtct["Option B", "Breastfeeding (per month) >=350"]
  ##This should be CD4 stratified, but only the <350 val is actually used
  mtct_trt[c('>500', '350-500'),"SDNVP","breastfeeding"] <- mtct["Single dose nevirapine", "Breastfeeding (per month) <350"]
  mtct_trt[c('250-349', '200-249', '100-199', '50-99', '<50'),"SDNVP","breastfeeding"] <- mtct["Single dose nevirapine", "Breastfeeding (per month) >=350"]
  mtct_trt[,"Dual ARV","breastfeeding"] <- mtct["WHO 2006 dual ARV regimen", "Breastfeeding (per month) >=350"]
  mtct_trt[c('250-349', '200-249', '100-199', '50-99', '<50'),"Option B+: before pregnancy","breastfeeding"] <- mtct["ART: Started before pregnancy", "Breastfeeding (per month) <350"]
  mtct_trt[c('250-349', '200-249', '100-199', '50-99', '<50'),"Option B+: >4 weeks","breastfeeding"] <- mtct["ART: Started during pregnancy >4 weeks", "Breastfeeding (per month) <350"]
  mtct_trt[c('250-349', '200-249', '100-199', '50-99', '<50'),"Option B+: <4 weeks","breastfeeding"] <- mtct["ART: Started during pregnancy <4 weeks", "Breastfeeding (per month) <350"]
  pmtct_mtct <- mtct_trt

  mtct <- array(data = 0, dim = c(8,2), dimnames = list(cd4 = c('>500', '350-500',
                                                                '250-349', '200-249',
                                                                '100-199', '50-99',
                                                                '<50', 'INFECTION'),
                                                        trans_type = c('perinatal', 'bf')))
  mtct[c('100-199', '50-99', '<50'),"perinatal"] <- untrt_mtct["CD4 <200","Perinatal"]
  mtct[c('250-349', '200-249'),"perinatal"] <- untrt_mtct["CD4 200-350","Perinatal"]
  mtct[c('>500', '350-500'),"perinatal"] <- untrt_mtct["CD4 >350","Perinatal"]
  mtct['INFECTION',"perinatal"] <- untrt_mtct["Incident infections","Perinatal"]

  mtct[c('100-199', '50-99', '<50'),"bf"] <- untrt_mtct["CD4 <200","Breastfeeding (per month) <350"]
  mtct[c('250-349', '200-249'),"bf"] <- untrt_mtct["CD4 200-350","Breastfeeding (per month) <350"]
  mtct[c('>500', '350-500'),"bf"] <- untrt_mtct["CD4 >350","Breastfeeding (per month) >=350"]
  mtct['INFECTION',"bf"] <- untrt_mtct["Incident infections","Breastfeeding (per month) >=350"]

  return(list(pmtct_mtct = pmtct_mtct, mtct = mtct))
}

prepare_hc_art_mort_rr <- function(data){
  ## projection parameters
  yr_start <- data$first_year$data
  yr_end <- data$final_year$data
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  year.idx <- 1:length(proj.years)

  mort_rr_art <- data$child_mortality_rates$data
  mort_rr_art_target <- array(NA, dim = c(3, 15, length(year.idx)),
                              dimnames = list(time_art = c('0to6mo', '7to12mo', '12+mo'),
                                              age = 0:14, year = proj.years))
  mort_rr_art_target[c('0to6mo', '7to12mo'),as.character(0:4),] <- rep(mort_rr_art["Age <5, <12 months on ART",],
                                                                       each = 10)
  mort_rr_art_target[c('12+mo'),as.character(0:4),] <- rep(mort_rr_art["Age <5, 12+ months on ART",],
                                                           each = 5)
  mort_rr_art_target[c('0to6mo', '7to12mo'),as.character(5:14),] <- rep(mort_rr_art["Age >=5, <12 months on ART",],
                                                                        each = 20)
  mort_rr_art_target[c('12+mo'),as.character(5:14),] <- rep(mort_rr_art["Age >=5, 12+ months on ART",],
                                                            each = 10)
  return(mort_rr_art_target)
}

prepare_art_elig <- function(data){
  ## projection parameters
  yr_start <- data$first_year$data
  yr_end <- data$final_year$data
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  year.idx <- 1:length(proj.years)

  hc_art_elig_age <- as.integer(data$age_hiv_child_on_treatment$data / 12) ##converts from months to years

  cd4_elig <- data$cd4_threshold$data[c("CD4 percent: <11 MOS",
                                        "CD4 percent: 12-35 MOS",
                                        "CD4 percent: 35-39 MOS",
                                        "CD4 count: >=5 YRS"),]
  ##Changing the input from CD4 count or percentages to ordinal categories
  ###Easier to do it here than in the leapfrog code
  paed_cd4_percent_intervals <- c(31, 30, 25, 20, 15, 10, 5)
  paed_cd4_number_intervals <- c(1001, 1000, 750, 500, 350, 200)
  hc_art_elig_cd4 <- array(data = NA, dim = c(15, length(year.idx)),
                             dimnames = list(age = c(0:14),
                                             year = c(proj.years)))
  hc_art_elig_cd4["0",] <- findInterval(-unname(cd4_elig["CD4 percent: <11 MOS",]), -paed_cd4_percent_intervals)
  hc_art_elig_cd4[as.character(1:2),] <- findInterval(-unname(cd4_elig["CD4 percent: 12-35 MOS",]), -paed_cd4_percent_intervals)
  hc_art_elig_cd4[as.character(3:4),] <- findInterval(-unname(cd4_elig["CD4 percent: 35-39 MOS",]), -paed_cd4_percent_intervals)
  hc_art_elig_cd4[as.character(5:14),] <- rep(findInterval(-unname(cd4_elig["CD4 count: >=5 YRS",]), -paed_cd4_number_intervals), each = length(6:15))

  return(list(hc_art_elig_age = hc_art_elig_age, hc_art_elig_cd4 = hc_art_elig_cd4))
}

prepare_bypass_adult_model <- function(data, bypass_adult = F){
  yr_start <- data$first_year$data
  yr_end <- data$final_year$data
  year.idx <- 1:length(yr_start:yr_end)

  ##extract needed outputs to just run paed model
  mat_hiv_births <- as.array(data$child_need_pmtct$data)
  ##Set this so that the default is to not use the input of maternal inputs, and instead use results from the adult model
  mat_prev_input = rep(bypass_adult, length(year.idx))

  hivpop <- data$hiv_by_single_age$data[as.character(15:49),"female",]
  totpop <- data$big_pop$data[as.character(15:49),"female",]
  inc <- data$new_infections_by_single_age$data[paste0(15:49, ' female'),]
  wlhiv_cd4 <- data$cd4_distribution_15_49$data[paste0("HIV, female: ", c("500+", "350-500", "250-349", "200-249", "100-199", "50-99", "50-")),]

  hivnpop <- totpop - hivpop
  adult_female_infections_full <- inc #/ hivnpop

  #cd4
  prop_gte350 <- colSums(wlhiv_cd4[paste0("HIV, female: ", c("500+", "350-500")),]) / colSums(wlhiv_cd4)
  prop_lt200 <- colSums(wlhiv_cd4[paste0("HIV, female: ", c("100-199", "50-99", "50-")),]) / colSums(wlhiv_cd4)

  return(list(mat_hiv_births = mat_hiv_births,
              mat_prev_input = mat_prev_input,
              hivnpop = hivnpop,
              adult_female_infections = adult_female_infections_full,
              prop_gte350 = prop_gte350,
              prop_lt200 = prop_lt200))
}

prepare_coarse_stratification <- function(data, params){
  adult_model <- prepare_bypass_adult_model(data)
  inc <- adult_model$adult_female_infections

  h.fert.idx <- which((15L-1 + cumsum(params$hAG_SPAN_coarse)) %in% 15:49)
  fert_rat.h.ag <- findInterval(as.integer(rownames(adult_model$hivnpop)),
                                15L + cumsum(params$hAG_SPAN_coarse[h.fert.idx]) -
                                  params$hAG_SPAN_coarse[h.fert.idx])
  coarse_age_groups <- cut(15:49,
                           breaks = c(15L + cumsum(params$hAG_SPAN_coarse[h.fert.idx]) -
                                        params$hAG_SPAN_coarse[h.fert.idx], 50),
                           right = FALSE)

  hivnpop_coarse <- as.array(rowsum(adult_model$hivnpop, group = coarse_age_groups))
  adult_female_infections_coarse <- as.array(rowsum(inc, group = coarse_age_groups))

  ##Make a coarse ASFR, needed to run WLHIV births at the coarse level
  fert_ages.idx <- 1L + cumsum(params$hAG_SPAN_coarse[h.fert.idx]) - params$hAG_SPAN_coarse[h.fert.idx]
  asfr_coarse <- as.array(rowsum(params$asfr, group = coarse_age_groups))

  return(list(hivnpop = hivnpop_coarse,
              adult_female_infections = adult_female_infections_coarse,
              asfr = asfr_coarse,
              fert_rat = params$fert_rat_coarse,
              frr_art6mos = params$frr_art6mos_coarse))
}

prepare_full_stratification <- function(data, params){
  adult_model <- prepare_bypass_adult_model(data)

  return(list(hivnpop = adult_model$hivnpop,
              adult_female_infections = adult_model$adult_female_infections,
              asfr = params$asfr,
              fert_rat = params$fert_rat_full,
              frr_art6mos = params$frr_art6mos_full))
}

#' Prepare child HIV projection parameters from Spectrum PJNZ
#'
#' @param pjnz path to PJNZ file
#' @param params The adult HIV projection parameters from [prepare_leapfrog_projp()]
#'
#' @return List of child HIV projection parameters
#'
#' @examples
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "frogger", mustWork = TRUE)
#' demp <- prepare_leapfrog_demp(pjnz)
#' proj <- prepare_leapfrog_projp(pjnz)
#' params <- c(proj, demp)
#' projp <- prepare_hc_leapfrog_projp(pjnz, params)
#' @export
prepare_hc_leapfrog_projp <- function(pjnz, params,
                                      bypass_adult = FALSE,
                                      use_coarse_age_groups = FALSE) {
  dp <- read_dp(pjnz)
  dat <- parse_dp(dp)
  data <- dat$data
  v <- params

  ## projection parameters
  yr_start <- data$first_year$data
  yr_end <- data$final_year$data
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  year.idx <- 1:length(proj.years)

  #############################################################
  ##Bypass adult inputs
  #############################################################
  bypass_adult <- prepare_bypass_adult_model(data, bypass_adult = bypass_adult)
  v$total_births <- data$births$data
  v$mat_hiv_births <- bypass_adult$mat_hiv_births
  v$mat_prev_input <- as.integer(bypass_adult$mat_prev_input)
  v$prop_gte350 <- bypass_adult$prop_gte350
  v$prop_lt200 <- bypass_adult$prop_lt200
  ###Need sex ratio of 1 and 2 year olds for option when running model without adult input
  v$infant_pop <- data$big_pop$data[as.character(1:2),,]

  #############################################################
  ##Births to WLHIV
  #############################################################
  if(use_coarse_age_groups){
    subparms = prepare_coarse_stratification(data, params = v)
  }else{
    subparms = prepare_full_stratification(data, params = v)
  }
  v$adult_female_hivnpop <- subparms$hivnpop
  v$adult_female_infections <- subparms$adult_female_infections
  v$hc_age_specific_fertility_rate <- subparms$asfr
  v$fert_mult_by_age <- subparms$fert_rat
  v$fert_mult_on_art <- subparms$frr_art6mos

  v$abortion <- prepare_abortion_input(data)
  v$patients_reallocated <- data$dp_tgx_patients_reallocated$data

  #############################################################
  ##Paediatric incidence
  #############################################################
  ##PMTCT
  pmtct <- prepare_pmtct(data)
  v$PMTCT <- pmtct$pmtct_new
  v$PMTCT_input_is_percent <- as.integer(pmtct$pmtct_input_isperc)
  v$PMTCT_dropout <- prepare_pmtct_dropout(data)

  ##rates of MTCT
  mtct <- prepare_vertical_transmission(data)
  v$PMTCT_transmission_rate <- mtct$pmtct_mtct
  v$vertical_transmission_rate <- mtct$mtct

  ##BF duration
  v$breastfeeding_duration_art <- data$infant_feeding_options$data[!grepl('No ART', rownames(data$infant_feeding_options$data)),] / 100
  v$breastfeeding_duration_no_art <- data$infant_feeding_options$data[grepl('No ART', rownames(data$infant_feeding_options$data)),] / 100

  ##TODO: change this in leapfrog to be by multiple ages
  v$hc_nosocomial <- data$nosocomial_infections_by_age$data[1,]

  v$hc1_cd4_dist <- data$child_dist_new_infections_cd4$data / 100

  #############################################################
  ##Natural history
  #############################################################
  prog <- prepare_cd4_progression(data)
  v$hc1_cd4_prog <- prog$hc1_cd4_prog
  v$hc2_cd4_prog <- prog$hc2_cd4_prog

  mort <- prepare_no_art_mort(data)
  v$hc1_cd4_mort <- mort$hc1_cd4_mort
  v$hc2_cd4_mort <- mort$hc2_cd4_mort

  #############################################################
  ##Paediatric treatment
  #############################################################
  ## Treatment eligibility
  art_elig = prepare_art_elig(data)
  v$hc_art_elig_age <- art_elig$hc_art_elig_age
  v$hc_art_elig_cd4 <- art_elig$hc_art_elig_cd4
  v$hc_art_init_dist <- data$child_art_dist$data

  ## Cotrim coverage
  v$ctx_val_is_percent <- ifelse(data$child_art_by_age_group_pernum$data["Cotrim",] == 1, T, F)
  v$ctx_val_is_percent <- as.integer(v$ctx_val_is_percent)
  v$ctx_val <- data$child_treat_inputs$data["Cotrim",]
  if(any(v$ctx_val_is_percent == 1)){
    v$ctx_val[v$ctx_val_is_percent == 1] <- v$ctx_val[v$ctx_val_is_percent == 1] / 100
  }
  v$ctx_effect <- prepare_cotrim_effect(data)

  ## ART coverage
  art = data$child_treat_inputs$data
  v$hc_art_is_age_spec <- unname(art[c("ART: 0-14y"),] == -9999)
  art[which(art == -9999)] <- 0
  ## only 0-14 input can be put in as percent
  v$hc_art_isperc <- as.integer(data$child_art_by_age_group_pernum$data["ART: 0-14y",])
  art["ART: 0-14y",which(v$hc_art_isperc == 1)] <- art["ART: 0-14y",which(v$hc_art_isperc == 1)] / 100
  v$hc_art_val <- art[c("ART: 0-14y", "ART: 0-4y", "ART: 5-9y", "ART: 10-14y"),]
  v$hc_art_start <- as.integer(unname(which(colSums(v$hc_art_val) > 0)[1]) - 1)

  v$hc_art_ltfu <- data$perc_lost_follow_up_child$data / 100

  #############################################################
  ##On ART mortality
  #############################################################
  art_mort <- prepare_art_mort(data)
  v$hc1_art_mort <- art_mort$hc1_art_mort
  v$hc2_art_mort <- art_mort$hc2_art_mort
  v$hc_art_mort_rr <- prepare_hc_art_mort_rr(data)


  return(v)
}

