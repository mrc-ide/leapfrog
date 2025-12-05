prepare_bypass_adult_model <- function(dat, pars, dim_vars, year_idx, bypass_adult) {
  ages_1549 <- as.character(15:49)
  cd4_greater_350 = c("500+", "350-500")
  cd4_less_200 = c("100-199", "50-99", "50-")

  ##extract needed outputs to just run paed model
  mat_hiv_births <- as.array(pars$child_need_pmtct)
  ##Set this so that the default is to not use the input of maternal inputs, and instead use results from the adult model
  mat_prev_input <- rep(bypass_adult, length(year_idx))

  hivpop <- pars$hiv_by_single_age[ages_1549, "female", ]
  totpop <- pars$big_pop[ages_1549, "female", ]
  inc <- pars$new_infections_by_single_age["female", ages_1549 , ]
  if (dat$data$cd4_distribution_15_49$tag == "CD4Distribution15_49 MV2") {
    wlhiv_cd4 <- pars$cd4_distribution_15_49[dim_vars$cd4_count, "female", "no art", ]
    prop_gte350 <- colSums(wlhiv_cd4[cd4_greater_350, ]) / colSums(wlhiv_cd4)
    prop_lt200 <-  colSums(wlhiv_cd4[cd4_less_200, ]) / colSums(wlhiv_cd4)
  } else {
    wlhiv_cd4 <- pars$cd4_distribution_15_49[paste0("HIV, female: ", dim_vars$cd4_count), ]
    prop_gte350 <- colSums(wlhiv_cd4[paste0("HIV, female: ", cd4_greater_350), ]) / colSums(wlhiv_cd4)
    prop_lt200 <-   colSums(wlhiv_cd4[paste0("HIV, female: ", cd4_less_200), ]) / colSums(wlhiv_cd4)
  }

  hivnpop <- totpop - hivpop
  adult_female_infections_full <- inc

  list(mat_hiv_births = mat_hiv_births,
       mat_prev_input = mat_prev_input,
       hivnpop = hivnpop,
       adult_female_infections = adult_female_infections_full,
       prop_gte350 = prop_gte350,
       prop_lt200 = prop_lt200)
}

prepare_coarse_stratification <- function(dat, pars, dim_vars) {
  inc <- pars$adult_female_infections

  h_fert_idx <- which((15L - 1 + cumsum(pars$h_ag_span_coarse)) %in% 15:49)
  fert_rat_h_ag <- findInterval(
    as.integer(rownames(pars$hivnpop)),
    15L + cumsum(pars$h_ag_span_coarse[h_fert_idx]) - pars$h_ag_span_coarse[h_fert_idx]
  )
  coarse_age_groups <- cut(
    15:49,
    breaks = c(
      15L + cumsum(pars$h_ag_span_coarse[h_fert_idx]) - pars$h_ag_span_coarse[h_fert_idx],
      50
    ),
    right = FALSE
  )

  hivnpop_coarse <- as.array(rowsum(pars$hivnpop, group = coarse_age_groups))
  adult_female_infections_coarse <- as.array(rowsum(inc, group = coarse_age_groups))

  ##Make a coarse ASFR, needed to run WLHIV births at the coarse level
  fert_ages_idx <- 1L + cumsum(pars$h_ag_span_coarse[h_fert_idx]) - pars$h_ag_span_coarse[h_fert_idx]
  asfr_coarse <- as.array(rowsum(pars$asfr, group = coarse_age_groups))

  list(hivnpop = hivnpop_coarse,
       adult_female_infections = adult_female_infections_coarse,
       asfr = asfr_coarse,
       fert_rat = pars$fert_rat_coarse,
       frr_art6mos = pars$frr_art6mos_coarse)
}

prepare_full_stratification <- function(dat, pars, dim_vars) {
  list(hivnpop = pars$hivnpop,
       adult_female_infections = pars$adult_female_infections,
       asfr = pars$asfr,
       fert_rat = pars$fert_rat_full,
       frr_art6mos = pars$frr_art6mos_full)
}

prepare_abortion_input <- function(dat, pars, dim_vars, proj_years) {
  abortion <- array(
    0, dim = c(2, length(proj_years)), dimnames = list(value = c("Input", "Percent"), year = proj_years)
  )
  abortion["Input", ] <- pars$preg_term_abortion
  abortion["Percent", ] <- pars$preg_term_abortion_pernum
  abortion
}

prepare_pmtct <- function(dat, pars, dim_vars) {
  pmtct <- pars$arv_regimen
  pmtct[is.na(pmtct)] <- 0
  pmtct_number <- pmtct[grepl("- Number", rownames(pmtct)) & !grepl("postnatal", rownames(pmtct))
                        & !grepl("Postnatal", rownames(pmtct)) & !grepl("Total", rownames(pmtct)), ]
  pmtct_pct <- pmtct[grepl("- Percent", rownames(pmtct)) & !grepl("postnatal", rownames(pmtct))
                     & !grepl("Postnatal", rownames(pmtct)) & !grepl("Total", rownames(pmtct)), ]
  pmtct_pct <- pmtct_pct[-(which(rownames(pmtct_pct) == "No prophylaxis- Percent")), ]
  pmtct_input_isperc <- rep(1, ncol(pmtct_number))
  pmtct_input_isperc[colSums(pmtct_number) > 0] <- 0
  pmtct_options <- c("Option A", "Option B", "Single dose nevirapine",
             "WHO 2006 dual ARV regimen", "ART: Started before pregnancy",
             "ART: Started during pregnancy >4 weeks", "ART: Started during pregnancy <4 weeks")
  pmtct_options_length = length(pmtct_options)

  rownames(pmtct_number) <- gsub(pattern = "- Number", replacement = "", rownames(pmtct_number))
  rownames(pmtct_pct) <- gsub(pattern = "- Percent", replacement = "", rownames(pmtct_pct))
  ##pmtct_options in the expected order for leapfrog
  pmtct_number <- pmtct_number[match(pmtct_options, rownames(pmtct_number)), ]
  pmtct_pct <- pmtct_pct[match(pmtct_options, rownames(pmtct_pct)), ]

  pmtct_new <- array(0, dim = c(pmtct_options_length, 61), dimnames = list(
    pmtct = pmtct_options
  ))
  ## pick out which ones were inserted as numbers
  pmtct_new <- pmtct_number
  ## pick out which ones were inserted as percent
  pmtct_new[, which(pmtct_input_isperc == 1)] <- pmtct_pct[, which(pmtct_input_isperc == 1)]

  list(pmtct_new = pmtct_new,
       pmtct_input_isperc = pmtct_input_isperc)
}

prepare_pmtct_dropout <- function(dat, pars, dim_vars, proj_years) {
  pmtct_options <- c("Option A", "Option B", "Single dose nevirapine",
                     "WHO 2006 dual ARV regimen", "ART: Started before pregnancy",
                     "ART: Started during pregnancy >4 weeks", "ART: Started during pregnancy <4 weeks")
  pmtct_dropout_eligible <- c("Option A",
                              "Option B",
                              "ART: Started before pregnancy",
                              "ART: Started during pregnancy >4 weeks",
                              "ART: Started during pregnancy <4 weeks")

  pmtct_dropout <- array(
    0, dim = c(length(pmtct_options), 3, length(proj_years)),
    dimnames = list(pmtct = pmtct_options,
                    drop_out_by = c("Delivery", "<12MOS breastfeeding", ">12MOS breastfeeding"),
                    year = proj_years)
  )
  pmtct_dropout[, "Delivery", ] <- 100
  pmtct_dropout["ART: Started before pregnancy", "Delivery", ] <- pars$percent_art_delivery["Percent already on ART retained at delivery", ]
  pmtct_dropout["ART: Started during pregnancy >4 weeks", "Delivery", ] <- pars$percent_art_delivery["Percent starting ART retained at delivery", ]
  pmtct_dropout[pmtct_dropout_eligible, "<12MOS breastfeeding", ] <- rep(pars$arv_regimen[rownames(pars$arv_regimen) %in%
                                                                                              c("Monthly dropout breastfeeding: ART 0-12 months breastfeeding"), ], each = 5)
  pmtct_dropout[pmtct_dropout_eligible, ">12MOS breastfeeding", ] <- rep(pars$arv_regimen[rownames(pars$arv_regimen) %in%
                                                                                              c("Monthly dropout breastfeeding: ART 12+ months breastfeeding"), ], each = 5)
  pmtct_dropout[is.na(pmtct_dropout)] <- 0
  pmtct_dropout <- pmtct_dropout / 100
  pmtct_dropout
}

prepare_vertical_transmission <- function(dat, pars, dim_vars) {
  pmtct_options <- c("Option A", "Option B", "Single dose nevirapine",
                     "WHO 2006 dual ARV regimen", "ART: Started before pregnancy",
                     "ART: Started during pregnancy >4 weeks", "ART: Started during pregnancy <4 weeks")
  cd4_greater_350 = c("500+", "350-500")
  cd4_less_200 = c("100-199", "50-99", "50-")
  cd4_200_to_350 = c("250-349", "200-249")
  cd4_less_350 = c(cd4_200_to_350, cd4_less_200)
  trans_types = c("perinatal", "breastfeeding")

  mtct <- pars$trans_eff_assump / 100
  untrt_mtct <- mtct[!rownames(mtct) %in% pmtct_options, ]
  mtct <- mtct[match(pmtct_options, rownames(mtct)), ]
  mtct_trt <- array(data = 0, dim = c(length(dim_vars$cd4_count), length(pmtct_options), length(trans_types)), dimnames = list(cd4 = dim_vars$cd4_count,
                                                                pmtct_reg = pmtct_options,
                                                                transmission_type = trans_types))
  mtct_trt[, "Option A", "perinatal"] <- mtct["Option A", "Perinatal"]
  mtct_trt[, "Option B", "perinatal"] <- mtct["Option B", "Perinatal"]
  mtct_trt[, "Single dose nevirapine", "perinatal"] <- mtct["Single dose nevirapine", "Perinatal"]
  mtct_trt[, "WHO 2006 dual ARV regimen", "perinatal"] <- mtct["WHO 2006 dual ARV regimen", "Perinatal"]
  mtct_trt[, "ART: Started before pregnancy", "perinatal"] <- mtct["ART: Started before pregnancy", "Perinatal"]
  mtct_trt[, "ART: Started during pregnancy >4 weeks", "perinatal"] <- mtct["ART: Started during pregnancy >4 weeks", "Perinatal"]
  mtct_trt[, "ART: Started during pregnancy <4 weeks", "perinatal"] <- mtct["ART: Started during pregnancy <4 weeks", "Perinatal"]

  mtct_trt[cd4_less_350, "Option A", "breastfeeding"] <- mtct["Option A", "Breastfeeding (per month) >=350"]
  mtct_trt[cd4_less_350, "Option B", "breastfeeding"] <- mtct["Option B", "Breastfeeding (per month) >=350"]
  ##This should be CD4 stratified, but only the <350 val is actually used
  mtct_trt[cd4_greater_350, "Single dose nevirapine", "breastfeeding"] <- mtct["Single dose nevirapine", "Breastfeeding (per month) <350"]
  mtct_trt[cd4_less_350, "Single dose nevirapine", "breastfeeding"] <- mtct["Single dose nevirapine", "Breastfeeding (per month) >=350"]
  mtct_trt[, "WHO 2006 dual ARV regimen", "breastfeeding"] <- mtct["WHO 2006 dual ARV regimen", "Breastfeeding (per month) >=350"]
  mtct_trt[cd4_less_350, "ART: Started before pregnancy", "breastfeeding"] <- mtct["ART: Started before pregnancy", "Breastfeeding (per month) <350"]
  mtct_trt[cd4_less_350, "ART: Started during pregnancy >4 weeks", "breastfeeding"] <- mtct["ART: Started during pregnancy >4 weeks", "Breastfeeding (per month) <350"]
  mtct_trt[cd4_less_350, "ART: Started during pregnancy <4 weeks", "breastfeeding"] <- mtct["ART: Started during pregnancy <4 weeks", "Breastfeeding (per month) <350"]
  pmtct_mtct <- mtct_trt

  mtct <- array(data = 0, dim = c(length(dim_vars$cd4_count) + 1, length(trans_types)), dimnames = list(cd4 = c(dim_vars$cd4_count, "INFECTION"),
                                                         trans_type = trans_types))
  mtct[cd4_less_200, "perinatal"] <- untrt_mtct["CD4 <200", "Perinatal"]
  mtct[cd4_200_to_350, "perinatal"] <- untrt_mtct["CD4 200-350", "Perinatal"]
  mtct[cd4_greater_350, "perinatal"] <- untrt_mtct["CD4 >350", "Perinatal"]
  mtct["INFECTION", "perinatal"] <- untrt_mtct["Incident infections", "Perinatal"]

  mtct[cd4_less_200, "breastfeeding"] <- untrt_mtct["CD4 <200", "Breastfeeding (per month) <350"]
  mtct[cd4_200_to_350, "breastfeeding"] <- untrt_mtct["CD4 200-350", "Breastfeeding (per month) <350"]
  mtct[cd4_greater_350, "breastfeeding"] <- untrt_mtct["CD4 >350", "Breastfeeding (per month) >=350"]
  mtct["INFECTION", "breastfeeding"] <- untrt_mtct["Incident infections", "Breastfeeding (per month) >=350"]

  list(pmtct_mtct = pmtct_mtct, mtct = mtct)
}

prepare_cd4_progression <- function(dat, pars, dim_vars) {
  hc2_cd4_cat = c(">1000", "750-999", "500-749", "350-499", "200-349", "lte200")

  prog <- pars$child_ann_rate_progress_lower_cd4
  hc1_cd4_prog <- array(0, dim = c(length(dim_vars$cd4_perc_0to4), 2, length(dim_vars$g)), dimnames = list(cd4pct = dim_vars$cd4_perc_0to4,
                                                             age = c("0-2", "3-4"),
                                                             sex = dim_vars$g))
  hc2_cd4_prog <- array(0, dim = c(length(hc2_cd4_cat), 1, length(dim_vars$g)), dimnames = list(cd4 = hc2_cd4_cat,
                                                             age = c("5-14"),
                                                             sex = dim_vars$g))
  hc1_cd4_prog[, "0-2", "male"] <- c(prog["male", grepl("Age: 0-2,", colnames(prog))], 0)
  hc1_cd4_prog[, "0-2", "female"] <- c(prog["female", grepl("Age: 0-2,", colnames(prog))], 0)
  hc1_cd4_prog[, "3-4", "male"] <- c(prog["male", grepl("Age: 3-4,", colnames(prog))], 0)
  hc1_cd4_prog[, "3-4", "female"] <- c(prog["female", grepl("Age: 3-4,", colnames(prog))], 0)
  hc2_cd4_prog[, "5-14", "male"] <- c(prog["male", grepl("Age: 5-14,", colnames(prog))], 0)
  hc2_cd4_prog[, "5-14", "female"] <- c(prog["female", grepl("Age: 5-14,", colnames(prog))], 0)

  list(hc1_cd4_prog = hc1_cd4_prog,
       hc2_cd4_prog = hc2_cd4_prog)
}

prepare_no_art_mort <- function(dat, pars, dim_vars) {
  hc2_cd4_cat = c(">1000", "750-999", "500-749", "350-499", "200-349", "lte200")
  ages_0to2 = as.character(0:2)
  ages_3to4 = as.character(3:4)
  ages_less5 = as.character(0:4)
  ages_5to14 = as.character(5:14)

  mort <- pars$child_mort_by_cd4_no_art
  hc1_cd4_mort <- array(data = 0, dim = c(length(dim_vars$cd4_perc_0to4), length(dim_vars$vt_time), length(ages_less5)), dimnames = list(cd4 = dim_vars$cd4_perc_0to4,
                                                                    transmission = dim_vars$vt_time,
                                                                    age = ages_less5))
  hc2_cd4_mort <- array(data = 0, dim = c(length(hc2_cd4_cat), length(dim_vars$vt_time), length(ages_5to14)), dimnames = list(cd4 = hc2_cd4_cat,
                                                                     transmission = dim_vars$vt_time,
                                                                     age = ages_5to14))
  ## 0-2
  hc1_cd4_mort[, "Perinatal", ages_0to2] <- mort["Perinatal", "00-02", ]
  hc1_cd4_mort[, "Breastfeeding, <6MOS", ages_0to2] <- mort["Breastfeeding, <6MOS", "00-02", ]
  hc1_cd4_mort[, "Breastfeeding, 7-12MOS", ages_0to2] <- mort["Breastfeeding, 7-12MOS", "00-02", ]
  hc1_cd4_mort[, "Breastfeeding, >12MOS", ages_0to2] <- mort["Breastfeeding, >12MOS", "00-02", ]

  ## 3-4
  hc1_cd4_mort[, "Perinatal", ages_3to4] <- mort["Perinatal", "03-04", ]
  hc1_cd4_mort[, "Breastfeeding, <6MOS", ages_3to4] <- mort["Breastfeeding, <6MOS", "03-04", ]
  hc1_cd4_mort[, "Breastfeeding, 7-12MOS", ages_3to4] <- mort["Breastfeeding, 7-12MOS", "03-04", ]
  hc1_cd4_mort[, "Breastfeeding, >12MOS", ages_3to4] <- mort["Breastfeeding, >12MOS", "03-04", ]

  ## 5-14, skip first index to account for the difference in CD4 categories
  hc2_cd4_mort[, "Perinatal", ages_5to14] <- mort["Perinatal", "05-14", ][-1]
  hc2_cd4_mort[, "Breastfeeding, <6MOS", ages_5to14] <- mort["Breastfeeding, <6MOS", "05-14", ][-1]
  hc2_cd4_mort[, "Breastfeeding, 7-12MOS", ages_5to14] <- mort["Breastfeeding, 7-12MOS", "05-14", ][-1]
  hc2_cd4_mort[, "Breastfeeding, >12MOS", ages_5to14] <- mort["Breastfeeding, >12MOS", "05-14", ] [-1]

  list(hc1_cd4_mort = hc1_cd4_mort,
       hc2_cd4_mort = hc2_cd4_mort)
}

prepare_art_elig <- function(dat, pars, dim_vars, proj_years, year_idx) {
  n_months = 12
  hc_art_elig_age <- as.integer(pars$age_hiv_child_on_treatment / n_months) ##converts from months to years

  cd4_elig <- pars$cd4_threshold[c("CD4 percent: <11 MOS",
                                   "CD4 percent: 12-35 MOS",
                                   "CD4 percent: 35-39 MOS",
                                   "CD4 count: >=5 YRS"), ]
  ##Changing the input from CD4 count or percentages to ordinal categories
  ###Easier to do it here than in the leapfrog code
  paed_cd4_percent_intervals <- c(31, 30, 25, 20, 15, 10, 5)
  paed_cd4_number_intervals <- c(1001, 1000, 750, 500, 350, 200)
  hc_art_elig_cd4 <- array(data = NA, dim = c(length(dim_vars$a_0to14), length(year_idx)),
                           dimnames = list(age = dim_vars$a_0to14,
                                           year = c(proj_years)))
  hc_art_elig_cd4["0", ] <- findInterval(-unname(cd4_elig["CD4 percent: <11 MOS", ]), -paed_cd4_percent_intervals)
  hc_art_elig_cd4[as.character(1:2), ] <- findInterval(-unname(cd4_elig["CD4 percent: 12-35 MOS", ]), -paed_cd4_percent_intervals)
  hc_art_elig_cd4[as.character(3:4), ] <- findInterval(-unname(cd4_elig["CD4 percent: 35-39 MOS", ]), -paed_cd4_percent_intervals)
  hc_art_elig_cd4[as.character(5:14), ] <- rep(findInterval(-unname(cd4_elig["CD4 count: >=5 YRS", ]), -paed_cd4_number_intervals), each = length(6:15))

  list(hc_art_elig_age = hc_art_elig_age, hc_art_elig_cd4 = hc_art_elig_cd4)
}

prepare_cotrim_effect <- function(dat, pars, dim_vars) {
  ##cotrim is effective for five years for children not on ART and for four years for children on ART
  ctx_yrs_no_art = 5
  ctx_yrs_art = 4

  ctx_effect <- pars$effect_treat_child
  off_art_ctx <- sum(as.numeric(unlist(ctx_effect["no art", ]))) / ctx_yrs_no_art
  on_art_ctx_lte12mo <- sum(as.numeric(unlist(ctx_effect["art", 1])))
  on_art_ctx_gte12mo <- sum(as.numeric(unlist(ctx_effect["art", 2:5]))) /  ctx_yrs_art
  ctx_effect <- array(data = c(off_art_ctx, on_art_ctx_lte12mo, on_art_ctx_gte12mo),
                      dim = c(3),
                      dimnames = list(ctx_effect = c("Off ART", "On ART, lte12mo", "On ART, gte12mo")))
  ctx_effect
}

prepare_art_mort <- function(dat, pars, dim_vars) {
  hc2_cd4_cat = c(">1000", "750-999", "500-749", "350-499", "200-349", "lte200")
  time_art = c("0to6mo", "7to12mo", "12+mo")
  ages_0to2 = as.character(0:2)
  ages_1to2 = as.character(1:2)
  ages_3to4 = as.character(3:4)
  ages_less5 = as.character(0:4)
  ages_5to9 = as.character(5:9)
  ages_10to14 = as.character(5:9)
  ages_5to14 = as.character(5:14)


  mort_lt6 <- pars$child_mort_by_cd4_with_art_0to6
  mort_7to12 <- pars$child_mort_by_cd4_with_art_7to12
  mort_gt12 <- pars$child_mort_by_cd4_with_art_gt12
  ##TODO: stratify this by sex
  hc1_art_mort <- array(data = 0, dim = c(length(dim_vars$cd4_perc_0to4), length(time_art), length(ages_less5)), dimnames = list(cd4 = dim_vars$cd4_perc_0to4,
                                                                    transmission = time_art,
                                                                    age = ages_less5))
  hc2_art_mort <- array(data = 0, dim = c(length(hc2_cd4_cat), length(time_art), length(ages_5to14)), dimnames = list(cd4 = hc2_cd4_cat,
                                                                     transmission = time_art,
                                                                     age = ages_5to14))
  ## 0-6 mo on treatment
  hc1_art_mort[, "0to6mo", "0"] <- mort_lt6["male", grepl("0: ", colnames(mort_lt6))]
  hc1_art_mort[, "0to6mo", ages_1to2] <- mort_lt6["male", grepl("1-2: ", colnames(mort_lt6))]
  hc1_art_mort[, "0to6mo", ages_3to4] <- mort_lt6["male", grepl("3-4: ", colnames(mort_lt6))]
  hc2_art_mort[, "0to6mo", ages_5to9] <- mort_lt6["male", grepl("5-9: ", colnames(mort_lt6))]
  hc2_art_mort[, "0to6mo", ages_10to14] <- mort_lt6["male", grepl("10-14: ", colnames(mort_lt6))]

  ## 7-12 mo on treatment
  hc1_art_mort[, "7to12mo", "0"] <- mort_7to12["male",grepl("0: ", colnames(mort_7to12))]
  hc1_art_mort[, "7to12mo", ages_1to2] <- mort_7to12["male", grepl("1-2: ", colnames(mort_7to12))]
  hc1_art_mort[, "7to12mo", ages_3to4] <- mort_7to12["male", grepl("3-4: ", colnames(mort_7to12))]
  hc2_art_mort[, "7to12mo", ages_5to9] <- mort_7to12["male", grepl("5-9: ", colnames(mort_7to12))]
  hc2_art_mort[, "7to12mo", ages_10to14] <- mort_7to12["male", grepl("10-14: ", colnames(mort_7to12))]

  ## 12+ mo on treatment
  hc1_art_mort[, "12+mo", "0"] <- mort_gt12["male", grepl("0: ", colnames(mort_gt12))]
  hc1_art_mort[, "12+mo", ages_1to2] <- mort_gt12["male", grepl("1-2: ", colnames(mort_gt12))]
  hc1_art_mort[, "12+mo", ages_3to4] <- mort_gt12["male", grepl("3-4: ", colnames(mort_gt12))]
  hc2_art_mort[, "12+mo", ages_5to9] <- mort_gt12["male", grepl("5-9: ", colnames(mort_gt12))]
  hc2_art_mort[, "12+mo", ages_10to14] <- mort_gt12["male", grepl("10-14: ", colnames(mort_gt12))]

  list(hc1_art_mort = hc1_art_mort,
       hc2_art_mort = hc2_art_mort)
}

prepare_hc_art_mort_rr <- function(dat, pars, dim_vars, proj_years, year_idx) {
  art_less_12mo = c("0to6mo", "7to12mo")
  art_more_12mo = "12+mo"
  time_art = c(art_less_12mo, art_more_12mo)
  ages_less5 = as.character(0:4)
  ages_5to14 = as.character(5:14)


  mort_rr_art <- pars$child_mortality_rates
  mort_rr_art_target <- array(NA, dim = c(length(time_art), length(dim_vars$a_0to14), length(year_idx)),
                              dimnames = list(time_art = time_art ,
                                              age = dim_vars$a_0to14, year = proj_years))
  mort_rr_art_target[art_less_12mo, ages_less5, ] <- rep(mort_rr_art["Age <5, <12 months on ART", ],
                                                                         each = 10)
  mort_rr_art_target[art_more_12mo, ages_less5, ] <- rep(mort_rr_art["Age <5, 12+ months on ART", ],
                                                             each = 5)
  mort_rr_art_target[art_less_12mo, ages_5to14, ] <- rep(mort_rr_art["Age >=5, <12 months on ART", ],
                                                                          each = 20)
  mort_rr_art_target[art_more_12mo, ages_5to14, ] <- rep(mort_rr_art["Age >=5, 12+ months on ART", ],
                                                              each = 10)
  mort_rr_art_target
}

process_pjnz_hc <- function(dat, pars, dim_vars, use_coarse_age_groups = FALSE, bypass_adult = FALSE) {
  ## projection parameters
  yr_start <- pars$first_year
  yr_end <- pars$final_year
  proj_years <- yr_start:yr_end
  timedat_idx <- 4 + seq_along(proj_years) - 1
  year_idx <- seq_along(proj_years)

  #############################################################
  # Bypass adult inputs
  #############################################################
  bypass_adult <- prepare_bypass_adult_model(dat, pars, dim_vars, year_idx, bypass_adult)
  pars$total_births <- pars$births
  pars$mat_hiv_births <- bypass_adult$mat_hiv_births
  pars$mat_prev_input <- as.integer(bypass_adult$mat_prev_input)
  pars$prop_gte350 <- bypass_adult$prop_gte350
  pars$prop_lt200 <- bypass_adult$prop_lt200
  ## Need sex ratio of 1 and 2 year olds for option when running model without adult input
  pars$infant_pop <- pars$big_pop[as.character(1:2), , ]
  pars$hivnpop <- bypass_adult$hivnpop
  pars$adult_female_infections <- bypass_adult$adult_female_infections

  #############################################################
  # Births to WLHIV
  #############################################################
  if (use_coarse_age_groups) {
    subparms <- prepare_coarse_stratification(dat, pars, dim_vars)
  } else {
    subparms <- prepare_full_stratification(dat, pars, dim_vars)
  }
  pars$adult_female_hivnpop <- subparms$hivnpop
  pars$adult_female_infections <- subparms$adult_female_infections
  pars$hc_age_specific_fertility_rate <- subparms$asfr
  pars$fert_mult_by_age <- subparms$fert_rat
  pars$fert_mult_on_art <- subparms$frr_art6mos

  pars$abortion <- prepare_abortion_input(dat, pars, dim_vars, proj_years)
  pars$patients_reallocated <- pars$dp_tgx_patients_reallocated

  #############################################################
  # Paediatric incidence
  #############################################################
  ##PMTCT
  pmtct <- prepare_pmtct(dat, pars, dim_vars)
  pars$PMTCT <- pmtct$pmtct_new
  pars$PMTCT_input_is_percent <- as.integer(pmtct$pmtct_input_isperc)
  pars$PMTCT_dropout <- prepare_pmtct_dropout(dat, pars, dim_vars, proj_years)

  ##rates of MTCT
  mtct <- prepare_vertical_transmission(dat, pars, dim_vars)
  pars$PMTCT_transmission_rate <- mtct$pmtct_mtct
  pars$vertical_transmission_rate <- mtct$mtct

  ##BF duration
  pars$breastfeeding_duration_art <- pars$infant_feeding_options[, "art", ] / 100
  pars$breastfeeding_duration_no_art <- pars$infant_feeding_options[, "no art", ] / 100

  ##TODO: change this in leapfrog to be by multiple ages
  pars$hc_nosocomial <- pars$nosocomial_infections_by_age[1, ]

  pars$hc1_cd4_dist <- pars$child_dist_new_infections_cd4 / 100

  #############################################################
  # Natural history
  #############################################################
  prog <- prepare_cd4_progression(dat, pars, dim_vars)
  pars$hc1_cd4_prog <- prog$hc1_cd4_prog
  pars$hc2_cd4_prog <- prog$hc2_cd4_prog

  mort <- prepare_no_art_mort(dat, pars, dim_vars)
  pars$hc1_cd4_mort <- mort$hc1_cd4_mort
  pars$hc2_cd4_mort <- mort$hc2_cd4_mort

  #############################################################
  # Paediatric treatment
  #############################################################
  ## Treatment eligibility
  art_elig <- prepare_art_elig(dat, pars, dim_vars, proj_years, year_idx)
  pars$hc_art_elig_age <- art_elig$hc_art_elig_age
  pars$hc_art_elig_cd4 <- art_elig$hc_art_elig_cd4
  pars$hc_art_init_dist <- pars$child_art_dist

  ## Cotrim coverage
  pars$ctx_val_is_percent <- pars$child_art_by_age_group_pernum["Cotrim", ] == 1
  pars$ctx_val_is_percent <- as.integer(pars$ctx_val_is_percent)
  pars$ctx_val <- pars$child_treat_inputs["Cotrim", ]
  if (any(pars$ctx_val_is_percent == 1)) {
    pars$ctx_val[pars$ctx_val_is_percent == 1] <- pars$ctx_val[pars$ctx_val_is_percent == 1] / 100
  }
  pars$ctx_effect <- prepare_cotrim_effect(dat, pars, dim_vars)

  ## ART coverage
  art <- pars$child_treat_inputs
  pars$hc_art_is_age_spec <- unname(art[c("ART: 0-14y"), ] == -9999)
  art[which(art == -9999)] <- 0
  ## only 0-14 input can be put in as percent
  pars$hc_art_isperc <- as.integer(pars$child_art_by_age_group_pernum["ART: 0-14y", ])
  art["ART: 0-14y", which(pars$hc_art_isperc == 1)] <- art["ART: 0-14y", which(pars$hc_art_isperc == 1)] / 100
  pars$hc_art_val <- art[c("ART: 0-14y", "ART: 0-4y", "ART: 5-9y", "ART: 10-14y"), ]
  pars$hc_art_start <- as.integer(unname(which(colSums(pars$hc_art_val) > 0)[1]) - 1)

  pars$hc_art_ltfu <- pars$perc_lost_follow_up_child / 100

  #############################################################
  # On ART mortality
  #############################################################
  art_mort <- prepare_art_mort(dat, pars, dim_vars)
  pars$hc1_art_mort <- art_mort$hc1_art_mort
  pars$hc2_art_mort <- art_mort$hc2_art_mort
  pars$hc_art_mort_rr <- prepare_hc_art_mort_rr(dat, pars, dim_vars, proj_years, year_idx)

  pars
}
