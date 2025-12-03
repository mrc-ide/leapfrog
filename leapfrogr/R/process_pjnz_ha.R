process_pjnz_ha <- function(dat, pars, dim_vars, dp, use_coarse_age_groups = FALSE) {
  if (!is.null(pars$incidence_by_fit) && !is.null(pars$incidence_options)) {
    pars$incidinput <- pars$incidence_by_fit[pars$incidence_options + 1, ]
  } else if (!is.null(pars$incidence_input)) {
    pars$incidinput <- pars$incidence_input
  } else {
    stop("Could not calculate incidinput")
  }
  pars$incidinput <- pars$incidinput / 100

  hiv_steps_per_year <- 10
  ds <- 7

  ## Incidence model inputs
  pars$incidence_model_choice <- 0L  ## 0: INCIDMOD_DIRECTINCID_HTS; 1: INCIDMOD_TRANSMRATE_HTS
  pars$transmission_rate_hts <- numeric(length(pars$incidinput) * hiv_steps_per_year)
  pars$initial_incidence <- 0.0
  pars$relative_infectiousness_art <- 0.1
  pars$epidemic_start_hts <- length(pars$transmission_rate_hts)

  pars$incrr_age <- pars$dist_of_hiv

  ## Note from Rob Glaubius (9 Dec 2020)
  ## Sometime between version 5.756 and the release v5.87 the variable <HIVSexRatio MV>
  ## was dropped for the more detailed <SexRatioByEpidPatt MV>. The old tag was
  ## reintroduced more recently with the shifts in cell placement you noticed. Since the
  ## variable was removed and reintroduced we didn’t think to tag this as “MV2” instead
  ## of “MV”.
  ##
  ## Action: to parse this check that one of these rows 2 or 3 after the <HIVSexRatio MV>
  ## has data, but not both.
  time_idx <- 4:(4 + length(dim_vars$years) - 1)
  if ("<HIVSexRatio MV>" %in% dp$Tag) {
    hiv_sex_ratio_idx <- which(dp[, 1] == "<HIVSexRatio MV>")
    first_entries <- dp[hiv_sex_ratio_idx + 2:3, 4]
    data_row <- which(first_entries != "" & !is.na(first_entries))
    if (length(data_row) != 1) {
      stop("DP tag <HIVSexRatio MV> not parsed correctly")
    }
    pars$incrr_sex <- array(
      as.numeric(unlist(dp[hiv_sex_ratio_idx + data_row + 1, time_idx])),
      dim = length(dim_vars$years),
      dimnames = list(dim_vars$years)
    )
  } else if ("<SexRatioByEpidPatt MV>" %in% dp$Tag) {
    sexincrr_idx <- as.integer(dp[which(dp[, 1] == "<IncEpidemicRGIdx MV>") + 2, 4]) # 0 based indexing
    data <- dp[which(dp[, 1] == "<SexRatioByEpidPatt MV>") + 3:8, time_idx][sexincrr_idx + 1, ]
    pars$incrr_sex <- array(
      as.numeric(data),
      dim = length(dim_vars$years),
      dimnames = list(dim_vars$years)
    )
  } else {
    stop("Incidence sex ratio not found")
  }

  ## Hard coded to expand age groups 15-24, 25-34, 35-44, 45+ to
  ## single-year ages 15:80.
  ## Requires extension for coarse HIV age group stratification
  idx_expand_full <- rep(1:4, times = c(10, 10, 10, 36))
  idx_expand_coarse <- rep(1:4, times = c(3, 2, 2, 2))
  idx_expand <- if (use_coarse_age_groups) idx_expand_coarse else idx_expand_full

  # reorder dims so sex is last dim
  pars$cd4_mort <- aperm(pars$adult_mort_by_cd4_no_art, c(2, 3, 1))
  pars$cd4_initdist <- aperm(pars$adult_dist_new_infections_cd4, c(2, 3, 1)) / 100
  pars$cd4_prog <- aperm(pars$adult_ann_rate_progress_lower_cd4, c(2, 3, 1))

  vers_str <- sub("^([0-9]+),(.*)$", "\\1.\\2", pars$valid_vers)
  version <- as.numeric(sub("^([0-9\\.]+).*", "\\1", vers_str))
  beta_version <- ifelse(
    grepl("Beta", vers_str),
    as.numeric(sub(".*Beta ([0-9]+)$", "\\1", vers_str)),
    NA
  )
  pars$scale_cd4_mort <- ifelse(
    version >= 5.73 && (beta_version >= 15 || is.na(beta_version)),
    1L,
    0L
  )

  ## eligibility starts in projection year idx
  pars$artcd4elig_idx <- findInterval(
    -pars$cd4_threshold_adults,
    -c(999, 500, 350, 250, 200, 100, 50)
  )
  # Update eligibility threshold from CD4 <200 to <250 to account for additional
  # proportion eligible with WHO Stage 3/4.
  pars$artcd4elig_idx <- replace(pars$artcd4elig_idx, pars$artcd4elig_idx == 5L, 4L)

  # create new dimension and combine arrays along it so now dims are:
  # sex, cd4 count, age groups, art duration
  pars$art_mort <- abind::abind(
    pars$adult_mort_by_cd4_with_art_0to6,
    pars$adult_mort_by_cd4_with_art_7to12,
    pars$adult_mort_by_cd4_with_art_gt12,
    along = 4
  )

  # dims are now: art duration, cd4 count, age groups, sex
  pars$art_mort <- aperm(pars$art_mort, c(4, 2, 3, 1))

  pars$artmx_timerr <- array(
    1, c(3, length(dim_vars$years)), list(
      artdur = c("ART0MOS", "ART6MOS", "ART1YR"),
      year = dim_vars$years
    )
  )
  if (!is.null(pars$mortality_rates)) {
    if (dat$data$mortality_rates$tag == "MortalityRates MV") {
      pars$artmx_timerr["ART0MOS", ] <- pars$mortality_rates
      pars$artmx_timerr["ART6MOS", ] <- pars$mortality_rates
      pars$artmx_timerr["ART1YR", ] <- pars$mortality_rates
    } else if (dat$data$mortality_rates$tag == "MortalityRates MV2") {
      pars$artmx_timerr["ART0MOS", ] <- pars$mortality_rates[1, ]
      pars$artmx_timerr["ART6MOS", ] <- pars$mortality_rates[1, ]
      pars$artmx_timerr["ART1YR", ] <- pars$mortality_rates[2, ]
    }
  }

  pars$cd4_nonaids_excess_mort <- array(0, dim(pars$cd4_mort), dimnames(pars$cd4_mort))
  pars$art_nonaids_excess_mort <- array(0, dim(pars$cd4_mort), dimnames(pars$cd4_mort))
  if (!is.null(pars$adult_non_aids_excess_mort)) {
    pars$cd4_nonaids_excess_mort <- aperm(pars$adult_non_aids_excess_mort[1, , , ], c(2, 3, 1))
    pars$art_nonaids_excess_mort <- aperm(pars$adult_non_aids_excess_mort[2, , , ], c(2, 3, 1))
  }

  pars$art_dropout_rate <- -log(1.0 - pars$perc_lost_followup / 100)

  pars$art15plus_numperc <- pars$ha_art_by_sex_per_num[c("male", "female"), ]
  pars$art15plus_num <- pars$ha_art_by_sex[c("male", "female"), ]

  adult_artadj_factor <- array(1, dim(pars$art15plus_num))
  adult_artadj_absolute <- array(0, dim(pars$art15plus_num))
  art_factor_flag <- is.null(pars$adult_art_adj_factor_flag) || pars$adult_art_adj_factor_flag == 1
  if (!is.null(pars$adult_art_adj_factor) && art_factor_flag) {
    adult_artadj_factor <- pars$adult_art_adj_factor
    if (!is.null(pars$adult_pats_alloc_to_from_other_region)) {
      adult_artadj_absolute <- pars$adult_pats_alloc_to_from_other_region
    }
    adult_artadj_factor <- adult_artadj_factor ^ as.numeric(!pars$art15plus_numperc)
    adult_artadj_absolute <- adult_artadj_absolute * as.numeric(!pars$art15plus_numperc)
    pars$art15plus_num <- pars$art15plus_num + adult_artadj_absolute
    pars$art15plus_num <- pars$art15plus_num * adult_artadj_factor
  }

  pars$art15plus_isperc <- pars$art15plus_numperc == 1
  pars$art15plus_num[pars$art15plus_isperc] <- pars$art15plus_num[pars$art15plus_isperc] / 100

  pars$art_alloc_mxweight <- pars$new_art_pat_alloc[1]

  p_ag_15to49 <- 35L
  p_ag_15plus <- 66L
  pars$pAG_INCIDPOP <- ifelse(pars$epp_population_ages == 0L, p_ag_15to49, p_ag_15plus)
  pars$pIDX_INCIDPOP <- 15L

  ## State space dimensions
  h_ag_span_full <- rep(1L, 66L)
  h_ag_span_coarse <- c(2L, 3L, 5L, 5L, 5L, 5L, 5L, 5L, 31L)
  pars$h_ag_span_coarse <- h_ag_span_coarse
  pars$h_ag_span_full <- h_ag_span_full

  ## Add in pediatric components
  frr_agecat <- seq(15, 45, 5)
  h_fert_idx <- which((15L - 1 + cumsum(h_ag_span_coarse)) %in% 15:49)
  age_band_width <- length(15:49) / length(frr_agecat)
  fert_rat_h_ag <- findInterval(
    15L + cumsum(h_ag_span_coarse[h_fert_idx]) - h_ag_span_coarse[h_fert_idx],
    frr_agecat
  )
  pars$fert_rat_coarse <- array(1, c(length(h_fert_idx), length(dim_vars$years)))
  pars$fert_rat_coarse[, ] <- pars$hivtfr[fert_rat_h_ag, ]
  pars$fert_rat_full <- apply(pars$hivtfr, 2, rep, each = age_band_width)

  if (is.null(pars$fert_cd4_discount)) {
    pars$cd4fert_rat <- rep(1, length(dim_vars$cd4_count))
  } else {
    pars$cd4fert_rat <- pars$fert_cd4_discount
  }

  if (is.null(pars$ratio_women_on_art)) {
    pars$ratio_women_on_art <- rep(1, length(dim_vars$a_15to49_5year))
  } else if (dat$data$ratio_women_on_art$tag == "RatioWomenOnART MV") {
    pars$ratio_women_on_art <- rep(pars$ratio_women_on_art, length(dim_vars$a_15to49_5year))
  }

  pars$frr_art6mos_full <- rep(pars$ratio_women_on_art, each = age_band_width)
  pars$frr_art6mos_coarse <- array(pars$ratio_women_on_art[fert_rat_h_ag])

  if (is.null(pars$frr_by_location)) {
    pars$frr_scalar <- 1
  } else {
    pars$frr_scalar <- pars$frr_by_location
  }

  if (use_coarse_age_groups) {
    pars$fert_rat <- pars$fert_rat_coarse
    pars$frr_art6mos <- pars$frr_art6mos_coarse
  } else {
    pars$fert_rat <- pars$fert_rat_full
    pars$frr_art6mos <- pars$frr_art6mos_full
  }

  proj_years <- as.integer(pars$final_year - pars$first_year + 1L)
  pars$t_ART_start <- min(c(unlist(apply(pars$art15plus_num > 0, 1, which)), proj_years))

  ## Use Beer's coefficients to distribution IRRs by age/sex
  pars$incrr_age <- apply(pars$incrr_age, 2:3, beers_open_ended)[16:81, , ] ## !! Hard coded
  pars$incrr_age[pars$incrr_age < 0] <- 0

  h_ts <- 3
  pars$cd4_initdist <- pars$cd4_initdist[, idx_expand, ]
  pars$cd4_prog <- (1 - exp(-pars$cd4_prog[, idx_expand, ] / hiv_steps_per_year)) * hiv_steps_per_year
  pars$cd4_mort <- pars$cd4_mort[, idx_expand, ]

  if (is.null(pars$mortality_rates_multiplier)) {
    pars$artmx_multiplier <- 1
  } else {
    pars$artmx_multiplier <- pars$mortality_rates_multiplier
  }
  pars$art_mort <- pars$artmx_multiplier * pars$art_mort[c(1, 2, 3), , idx_expand, ]
  pars$cd4_nonaids_excess_mort <- pars$cd4_nonaids_excess_mort[, idx_expand, ]
  art_nonaids_excess_mort_hts <- array(0.0, dim(pars$art_mort), dimnames(pars$art_mort))
  art_nonaids_excess_mort_hts[] <- rep(pars$art_nonaids_excess_mort[, idx_expand, ], each = h_ts)
  pars$art_nonaids_excess_mort <- art_nonaids_excess_mort_hts

  pars$art_dropout_recover_cd4 <- vers_str >= "6.14"

  pars
}
