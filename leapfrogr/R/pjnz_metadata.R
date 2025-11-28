get_data_from_cfg <- function(name, cfg, dim_vars, dp) {
  data <- NULL
  matching_tag_cfg <- NULL
  for (tag_cfg in cfg$read) {
    data <- get_data_from_tag_cfg(tag_cfg, dim_vars, dp)
    matching_tag_cfg <- tag_cfg
    if (!is.null(data)) break()
  }
  if (is.null(data)) {
    if (!is.null(cfg$allow_null) && cfg$allow_null) {
      warning(sprintf("Tag not found in DP for %s, returning NULL", name))
      return(NULL)
    } else {
      stop(sprintf("No tag recognised for %s", name))
    }
  }

  conversion_function <- function(x) x
  if (!is.null(cfg$type)) {
    if (cfg$type == "real") {
      conversion_function <- as.numeric
    } else if (cfg$type == "int") {
      conversion_function <- as.integer
    }
  }

  shaped_data <- NULL
  if (is.null(matching_tag_cfg$dims)) {
    shaped_data <- conversion_function(data)
  } else {
    dims <- dim_vars[unlist(matching_tag_cfg$dims)]
    dim_lengths <- unlist(lapply(dims, `[[`, "length"))
    dim_names <- lapply(dims, function(d) unlist(d$labels))
    shaped_data <- array(
      conversion_function(unlist(data)),
      dim = dim_lengths,
      dimnames = unname(dim_names)
    )
  }
  print(name)
  list(data = shaped_data, tag = matching_tag_cfg$tag)
}

get_data_from_tag_cfg <- function(tag_cfg, dim_vars, dp) {
  tag_idx <- which(dp[, 1] == sprintf("<%s>", tag_cfg$tag))
  if (length(tag_idx) == 0) return(NULL)

  # Not a fixed offset between different tags and version
  value_idx <- tag_idx
  while (dp$Description[value_idx] != "<Value>") value_idx <- value_idx + 1

  end_idx <- value_idx
  while (dp$Tag[end_idx] != "<End>") end_idx <- end_idx + 1

  start_row <- value_idx
  # we do not count the row with the <End> tag
  end_row <- end_idx - 1

  # PJNZ files have headers:
  # Tag, Description, Notes, Data, X, X.1, X.2, ...
  # we always want to start at Data
  start_column <- which(colnames(dp) == "Data")

  # some datasets have labels in first row such as years so we allow
  # metadata to specify custom start offset
  if (!is.null(tag_cfg$start_offset)) {
    offset <- tag_cfg$start_offset
    start_offset_row <- if (is.null(offset$row)) 0 else offset$row
    start_offset_column <- if (is.null(offset$column)) 0 else offset$column

    start_row <- start_row + start_offset_row
    start_column <- start_column + start_offset_column
  }

  # if no dims then we have a scalar so column end is 1
  # otherwise if user has specified the column_dims use the product of
  # those dim vars or assume the last dim is the length of the columns
  # needed to parse from the PJNZ
  n_cols <- if (is.null(tag_cfg$dims)) {
    1
  } else {
    if (!is.null(tag_cfg$column_dims)) {
      prod(unlist(lapply(tag_cfg$column_dims, function(var) dim_vars[[var]]$length)))
    } else {
      last_dim <- tag_cfg$dims[[length(tag_cfg$dims)]]
      dim_vars[[last_dim]]$length
    }
  }
  end_column <- start_column + n_cols - 1

  skip_rows <- NULL
  skip_columns <- NULL
  if (!is.null(tag_cfg$skip)) {
    if (!is.null(tag_cfg$skip$rows)) {
      # we do not add to the end row index because the end row is
      # calculated by wherever the <End> tag is which will include
      # the skipped rows
      skip_rows <- unlist(tag_cfg$skip$rows)
    }
    if (!is.null(tag_cfg$skip$columns)) {
      end_column <- end_column + length(tag_cfg$skip$columns)
      skip_columns <- unlist(tag_cfg$skip$columns)
    }
  }
  data_rows <- start_row:end_row
  data_columns <- start_column:end_column
  data_rectangle <- dp[data_rows, data_columns]

  if (!is.null(skip_rows)) {
    data_rectangle <- data_rectangle[-skip_rows, ]
  }
  if (!is.null(skip_columns)) {
    data_rectangle <- data_rectangle[, -skip_columns]
  }

  data_rectangle
}

get_dim_vars <- function(dp) {
  dim_vars <- get_static_dim_vars()
  add_years_to_dim_vars(dim_vars, dp)
}

add_years_to_dim_vars <- function(dim_vars, dp) {
  years_cfg <- get_years_cfg()
  start <- get_data_from_cfg("first_year", years_cfg$first_year, dim_vars, dp)$data
  end <- get_data_from_cfg("final_year", years_cfg$final_year, dim_vars, dp)$data

  dim_vars$years <- list(length = end - start + 1, labels = as.character(start:end))
  dim_vars
}

parse_dp <- function(dp) {
  dim_vars <- get_dim_vars(dp)
  metadata <- c(get_years_cfg(), get_pars_metadata(dim_vars, dp))

  ret <- lapply(names(metadata), function(name) {
    get_data_from_cfg(name, metadata[[name]], dim_vars, dp)
  })
  names(ret) <- names(metadata)
  list(data = ret, dim_vars = dim_vars)
}

get_pars_metadata <- function(dim_vars, dp) {
  list(
    version_num = list(
      read = list(
        list(tag = "VersionNum MV"),
        list(tag = "VersionNum MV2")
      )
    ),
    valid_vers = list(
      read = list(
        list(tag = "ValidVers MV")
      )
    ),
    big_pop = list(
      type = "real",
      read = list(
        list(
          tag = "BigPop MV",
          dims = list("a", "g", "years"),
          start_offset = list(row = 1)
        ),
        list(
          tag = "BigPop MV2",
          dims = list("a", "g_both_middle", "years"),
          start_offset = list(row = 1)
        ),
        list(
          tag = "BigPop MV3",
          dims = list("a", "g", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    surv_rate = list(
      type = "real",
      read = list(
        list(
          tag = "SurvRate MV",
          dims = list("a_80plus", "g", "years"),
          start_offset = list(row = 1),
          skip = list(
            rows = list(
              dim_vars$a_80plus$length + 1,
              dim_vars$a_80plus$length * 2 + 2
            )
          )
        ),
        list(
          tag = "SurvRate MV2",
          dims = list("a_80plus", "g", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    tfr = list(
      type = "real",
      read = list(
        list(
          tag = "TFR MV",
          dims = list("years")
        )
      )
    ),
    asfr = list(
      type = "real",
      read = list(
        list(
          tag = "ASFR MV",
          dims = list("a_15to49_5year", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    births = list(
      type = "real",
      read = list(
        list(
          tag = "Births MV",
          dims = list("years")
        )
      )
    ),
    sex_birth_ratio = list(
      type = "real",
      read = list(
        list(
          tag = "SexBirthRatio MV",
          dims = list("years")
        )
      )
    ),
    migr_rate = list(
      type = "real",
      read = list(
        list(
          tag = "MigrRate MV",
          dims = list("g", "years"),
          start_offset = list(row = 1),
          skip = list(
            # 3 and 6 have actual data, the rest are descriptions so
            # we have 2 bits of description above male and female
            rows = list(
              1, 2, 4, 5
            )
          )
        ),
        list(
          tag = "MigrRate MV2",
          dims = list("g", "years"),
          start_offset = list(row = 1),
          skip = list(
            # 2 and 4 have actual data,
            # we have 1 bit of description above male and female
            rows = list(1, 3)
          )
        )
      )
    ),
    migr_age_dist = list(
      type = "real",
      read = list(
        list(
          tag = "MigrAgeDist MV",
          dims = list("a_0to80plus_5year", "g", "years"),
          skip = list(
            rows = as.list(
              c(1:2, 3, 5,
                4 + 0:dim_vars$a_0to80plus_5year$length * 2,
                4 + dim_vars$a_0to80plus_5year$length * 2 + 2,
                4 + dim_vars$a_0to80plus_5year$length * 2 + 4,
                4 + dim_vars$a_0to80plus_5year$length * 2 + 2 + 0:dim_vars$a_0to80plus_5year$length * 2 + 1)
            )
          )
        ),
        list(
          tag = "MigrAgeDist MV2",
          dims = list("a_0to80plus_5year", "g", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    # DP DONE, START OF HA
    adult_infect_reduc = list(
      type = "real",
      read = list(
        list(
          tag = "AdultInfectReduc MV"
        )
      )
    ),
    hivtfr = list(
      type = "real",
      read = list(
        list(
          tag = "HIVTFR MV",
          dims = list("a_15to49_5year", "years")
        ),
        list(
          tag = "HIVTFR MV2",
          dims = list("a_15to35_custom", "years")
        ),
        list(
          tag = "HIVTFR MV3",
          dims = list("a_15to49_5year", "years")
        ),
        list(
          tag = "HIVTFR MV4",
          dims = list("a_15to49_5year", "years")
        )
      )
    ),
    fert_cd4_discount = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "FertCD4Discount MV",
          dims = list("cd4_count"),
          start_offset = list(column = 1)
        )
      )
    ),
    ratio_women_on_art = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "RatioWomenOnART MV"
        ),
        list(
          tag = "RatioWomenOnART MV2",
          dims = list("a_15to49_5year")
        )
      )
    ),
    frr_by_location = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "FRRbyLocation MV"
        )
      )
    ),
    dist_of_hiv = list(
      type = "real",
      read = list(
        list(
          tag = "DistOfHIV MV",
          dims = list("a_0to80plus_5year", "g", "years"),
          start_offset = list(row = 1),
          skip = list(
            # Males and Females label so ignore those rows
            rows = list(1, dim_vars$a_0to80plus_5year$length + 2)
          )
        ),
        list(
          tag = "DistOfHIV MV2",
          dims = list("a_0to80plus_5year", "g", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    adult_dist_new_infections_cd4 = list(
      type = "real",
      read = list(
        list(
          tag = "AdultDistNewInfectionsCD4 MV",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          start_offset = list(row = 1),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        )
      )
    ),
    adult_ann_rate_progress_lower_cd4 = list(
      type = "real",
      read = list(
        list(
          tag = "AdultAnnRateProgressLowerCD4 MV",
          dims = list("g", "cd4_prog", "a_15to45plus_10year"),
          start_offset = list(row = 1),
          column_dims = list("cd4_prog", "a_15to45plus_10year"),
          skip = list(
            columns = list(
              7, 14, 21, 28
            )
          )
        )
      )
    ),
    adult_mort_by_cd4_no_art = list(
      type = "real",
      read = list(
        list(
          tag = "AdultMortByCD4NoART MV",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          start_offset = list(row = 1),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        )
      )
    ),
    adult_mort_by_cd4_with_art_0to6 = list(
      type = "real",
      read = list(
        list(
          tag = "AdultMortByCD4WithART0to6 MV",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          start_offset = list(row = 1),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        ),
        list(
          tag = "AdultMortByCD4WithART0to6 MV2",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        )
      )
    ),
    adult_mort_by_cd4_with_art_7to12 = list(
      type = "real",
      read = list(
        list(
          tag = "AdultMortByCD4WithART7to12 MV",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          start_offset = list(row = 1),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        ),
        list(
          tag = "AdultMortByCD4WithART7to12 MV2",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        )
      )
    ),
    adult_mort_by_cd4_with_art_gt12 = list(
      type = "real",
      read = list(
        list(
          tag = "AdultMortByCD4WithARTGt12 MV",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          start_offset = list(row = 1),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        ),
        list(
          tag = "AdultMortByCD4WithARTGt12 MV2",
          dims = list("g", "cd4_count", "a_15to45plus_10year"),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        )
      )
    ),
    mortality_rates = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "MortalityRates MV",
          dims = list("years")
        ),
        list(
          tag = "MortalityRates MV2",
          dims = list("art_dur_2", "years")
        )
      )
    ),
    mortality_rates_multiplier = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "MortalityRatesMultiplier MV"
        )
      )
    ),
    adult_non_aids_excess_mort = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "AdultNonAIDSExcessMort MV",
          dims = list("art_status", "g", "cd4_count", "a_15to45plus_10year"),
          column_dims = list("cd4_count", "a_15to45plus_10year")
        )
      )
    ),
    ha_art_by_sex_per_num = list(
      type = "real",
      read = list(
        list(
          tag = "HAARTBySexPerNum MV",
          dims = list("g_both_first", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    ha_art_by_sex = list(
      type = "real",
      read = list(
        list(
          tag = "HAARTBySex MV",
          dims = list("g_both_first", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    cd4_threshold_adults = list(
      type = "real",
      read = list(
        list(
          tag = "CD4ThreshHoldAdults MV",
          dims = list("years")
        )
      )
    ),
    pops_elig_treat = list(
      type = "real",
      read = list(
        list(
          tag = "PopsEligTreat MV",
          dims = list("pops_for_treat", "elig_perc_year")
        )
      )
    ),
    valid_vers = list(
      read = list(
        list(
          tag = "ValidVers MV"
        )
      )
    ),
    med_cd4_count_init = list(
      type = "real",
      read = list(
        list(
          tag = "MedCD4CountInit MV",
          dims = list("years")
        )
      )
    ),
    perc_lost_followup = list(
      type = "real",
      read = list(
        list(
          tag = "PercLostFollowup MV",
          dims = list("years")
        )
      )
    ),
    hiv_by_single_age = list(
      type = "real",
      read = list(
        list(
          tag = "HIVBySingleAge MV",
          dims = list("a", "g", "years"),
          start_offset = list(row = 1),
          skip = list(
            rows = list(
              dim_vars$a$length + 1,
              dim_vars$a$length * 2 + 2
            )
          )
        ),
        list(
          tag = "HIVBySingleAge MV2",
          dims = list("a", "g", "years"),
          start_offset = list(row = 1)
        )
      )
    ),
    aids_deaths_by_age = list(
      type = "real",
      read = list(
        list(
          tag = "AidsDeathsByAge MV",
          dims = list("a", "g", "years"),
          start_offset = list(row = 2),
          skip = list(
            rows = list(
              dim_vars$a$length + 1
            )
          )
        ),
        list(
          tag = "AidsDeathsByAge MV2",
          dims = list("a", "g", "years"),
          start_offset = list(row = 1)
        )
      )
    ),

    # These three are all linked to incidence input
    incidence_by_fit = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "IncidenceByFit MV4",
          dims = list("incid_type", "years")
        )
      )
    ),
    incidence_options = list(
      allow_null = TRUE,
      type = "int",
      read = list(
        list(
          tag = "IncidenceOptions MV"
        )
      )
    ),
    incidence_input = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "IncidenceInput MV",
          dims = list("years")
        )
      )
    ),
    epp_population_ages = list(
      type = "int",
      read = list(
        list(
          tag = "EPPPopulationAges MV"
        )
      )
    ),

    adult_art_adj_factor = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "AdultARTAdjFactor",
          dims = list("years"),
          start_offset = list(row = 1)
        )
      )
    ),
    adult_art_adj_factor_flag = list(
      allow_null = TRUE,
      type = "int",
      read = list(
        list(
          tag = "AdultARTAdjFactorFlag"
        )
      )
    ),
    adult_pats_alloc_to_from_other_region = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "AdultPatsAllocToFromOtherRegion",
          dims = list("years"),
          start_offset = list(row = 1)
        )
      )
    ),
    new_art_pat_alloc = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "NewARTPatAlloc MV",
          dims = list("alloc_type")
        )
      )
    ),

    #########CHILD TAGS HEREAFTER
    trans_eff_assump = list(
      type = "real",
      read = list(
        list(
          tag = "TransEffAssump MV",
          dims = list("vt_trt", "vt_trt_cat")
        )
      )
    ),
    arv_regimen = list(
      type = "real",
      read = list(
        list(
          tag = "ARVRegimen MV3",
          dims = list("pmtct_editor_order", "years"),
          start_offset = list(column = 1)
        ),
        list(
          tag = "ARVRegimen MV2",
          dims = list("pmtct_editor_order_old", "years"),
          start_offset = list(column = 1)
        ),
        list(
          tag = "ARVRegimen MV",
          dims = list("pmtct_editor_order_old", "years"),
          start_offset = list(column = 1),
          skip = list(
            rows = as.list(
              c(1,
                3 + 0:6 * 3,
                24, 26:27,
                29 + 0:1 * 3,
                35, 37)
            )
          )
        )
      )
    ),
    percent_art_delivery = list(
      ##Note: not listed in spectrum files before 2018
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "PercentARTDelivery MV",
          dims = list("pmtct_retained_delivery", "years")
        )
      )
    ),
    preg_term_abortion = list(
      type = "real",
      read = list(
        list(
          tag = "PregTermAbortion MV3",
          dims = list("years")
        ),
        list(
          tag = "PregTermAbortion MV2",
          dims = list("years")
        ),
        list(
          tag = "PregTermAbortion MV",
          dims = list("years")
        )
      )
    ),
    preg_term_abortion_pernum = list(
      type = "real",
      read = list(
        list(
          tag = "PregTermAbortionPerNum MV2",
          dims = list("years")
        ),
        list(
          tag = "PregTermAbortionPerNum MV",
          dims = list("years")
        )
      )
    ),
    dp_tgx_patients_reallocated = list(
      ##Note: not listed in spectrum files before 2018
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "DP_TGX_PatientsReallocated_MV",
          dims = list("years")
        )
      )
    ),
    infant_feeding_options = list(
      type = "real",
      read = list(
        list(
          tag = "InfantFeedingOptions MV",
          dims = list("bf_months", "art_status", "years")
        )
      )
    ),
    child_treat_inputs = list(
      type = "real",
      read = list(
        list(
          tag = "ChildTreatInputs MV3",
          dims = list("paed_treatment", "years")
        ),
        list(
          tag = "ChildTreatInputs MV",
          dims = list("paed_treatment_old", "years")
        ),
        list(
          tag = "ChildTreatInputs2 MV",
          dims = list("paed_treatment", "years")
        )
      )
    ),
    child_art_by_age_group_pernum = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "ChildARTByAgeGroupPerNum MV2",
          dims = list("paed_treatment", "years")
        ),
        list(
          tag = "ChildARTByAgeGroupPerNum MV",
          dims = list("paed_treatment_old", "years")
        )
      )
    ),
    perc_lost_follow_up_child = list(
      ##Note: not listed in spectrum files before 2017
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "PercLostFollowupChild MV",
          dims = list("years")
        )
      )
    ),
    child_art_dist = list(
      type = "real",
      read = list(
        list(
          tag = "ChildARTDist MV",
          dims = list("a_0to14", "years")
        )
      )
    ),
    child_mortality_rates = list(
      ##Note: not listed in spectrum files before 2020
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "ChildMortalityRates MV2",
          dims = list("paed_art_mort", "years")
        ),
        list(
          tag = "ChildMortalityRates MV",
          dims = list("paed_art_mort", "years")
        )
      )
    ),
    ##TODO: Maggie is this implemented in leapfrog?
    child_mortality_rates_multiplier = list(
      ##Note: not listed in spectrum files before 2022
      allow_null = TRUE,
      type = "real",
      read = list(
        list(tag = "ChildMortalityRatesMultiplier MV")
      )
    ),
    nosocomial_infections_by_age = list(
      ##Note: not listed in spectrum files before 2020
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "NosocomialInfectionsByAge MV",
          dims = list("a_0to14_coarse", "years")
        )
      )
    ),
    nosocomial_infections = list(
      ##Note: not listed in spectrum files before 2018 or after 2020
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "NosocomialInfections MV",
          dims = list("years")
        )
      )
    ),
    child_dist_new_infections_cd4 = list(
      type = "real",
      read = list(
        list(
          tag = "ChildDistNewInfectionsCD4 MV",
          dims = list("cd4_perc_0to4")
        )
      )
    ),
    child_ann_rate_progress_lower_cd4 = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "ChildAnnRateProgressLowerCD4 MV2",
          dims = list("g", "cd4_prog_0to14"),
          skip = list(
            columns = (dim_vars$cd4_perc_0to4$length - 1) * 2 + 1
          ),
          start_offset = list(row = 1)
        ),
        list(
          tag = "ChildAnnRateProgressLowerCD4 MV",
          dims = list("g", "cd4_prog_0to14_old"),
          skip = list(
            column = (dim_vars$cd4_perc_0to4$length) * 2 + 1
          ),
          start_offset = list(row = 1)
        )
      )
    ),
    child_mort_by_cd4_no_art = list(
      allow_null = TRUE,
      type = "real",
      read = list(
        list(
          tag = "ChildMortByCD4NoART MV2",
          dims = list("vt_time", "a_0to14_coarse_2", "paed_ordinal_cd4_cat")
        ),
        ##Note: NOT DONE
        list(
          tag = "ChildMortByCD4NoART MV",
          dims = list("vt_time_old", "a_0to14_coarse_2", "paed_ordinal_cd4_cat", "vt_time_old"),
          column_dims = list("paed_ordinal_cd4_cat", "vt_time_old"),
          skip = list(
            rows = c(1,5:6,10:11,15:16)
          )
        )
      )
    ),
  child_mort_by_cd4_with_art_0to6 = list(
    type = "real",
    read = list(
      list(
        tag = "ChildMortByCD4WithART0to6 MV2",
        dims = list("g", "paed_age_time_art_cd4"),
        skip = list(
          columns = c((dim_vars$cd4_perc_0to4$length * 3) + 1,
                     (dim_vars$cd4_perc_0to4$length * 4) + 1)
        )
      ),
      list(
        tag = "ChildMortByCD4WithART0to6 MV",
        dims = list("g", "paed_age_time_art_cd4"),
        start_offset = list(row = 1),
        skip = list(
          columns = c((dim_vars$cd4_perc_0to4$length * 3) + 1,
                      (dim_vars$cd4_perc_0to4$length * 4) + 1)
        )
      )
    )
  ),
  child_mort_by_cd4_with_art_7to12 = list(
    type = "real",
    read = list(
      list(
        tag = "ChildMortByCD4WithART7to12 MV",
        dims = list("g", "paed_age_time_art_cd4"),
        start_offset = list(row = 1),
        skip = list(
          columns = c((dim_vars$cd4_perc_0to4$length * 3) + 1,
                      (dim_vars$cd4_perc_0to4$length * 4) + 1)
        )
      )
    )
  ),
  child_mort_by_cd4_with_art_gt12 = list(
    type = "real",
    read = list(
      list(
        tag = "ChildMortByCD4WithARTGT12 MV",
        dims = list("g", "paed_age_time_art_cd4"),
        start_offset = list(row = 1),
        skip = list(
          columns = c((dim_vars$cd4_perc_0to4$length * 3) + 1,
                      (dim_vars$cd4_perc_0to4$length * 4) + 1)
        )
      )
    )
  ),
 cd4_threshold = list(
   type = "real",
   read = list(
     list(
       tag = "CD4ThreshHold MV",
       dims = list("paed_cd4_art_elig", "years")
     )
   )
 ),
 age_hiv_child_on_treatment = list(
   type = "real",
   read = list(
     list(
       tag = "AgeHIVChildOnTreatment MV",
       dims = list("years")
     )
   )
 ),
 effect_treat_child = list(
   type = "real",
   read = list(
     list(
       tag = "EffectTreatChild MV",
       dims = list("art_status", "cotrim_years_effective"),
       start_offset = list(row = 1)
     )
   )
 ),
 child_need_pmtct = list(
   type = "real",
   read = list(
     list(
       tag = "ChildNeedPMTCT MV",
       dims = list("years")
     )
   )
 ),
##ASK: This is 15-49 year olds, not sure if I should make this clear in the labels
 cd4_distribution = list(
   allow_null = TRUE,
   type = "real",
   read = list(
     list(
       tag = "CD4Distribution MV2",
       dims = list("neg_adult_cd4_categories", "g_both_first", "art_status", "years"),
       skip = list(
         rows = c(25,50)
       )
       ##"CD4Distribution MV" exists, but the interpretation is difficult,
       ##seems like the first year is the year before HIV is introduced?
       ##could also be that HIV starts earlier in this file? but I cannot open this file in spectrum
       ##"C:/Users/mwalters/Imperial College London/HIV Inference Group - WP - Documents/Data/Spectrum files/2016 final/subnational/Malawi 2016eBlantyre.PJNZ"
     )
   )
 ),
new_infections_by_single_age = list(
  ##Not available before 2018
  allow_null = TRUE,
  type = "real",
  read = list(
    list(
      tag = "NewInfectionsBySingleAge MV",
      dims = list("g_both_first", "a_80plus_no_80", "years")
    )
  )
),
cd4_distribution_15_49 = list(
  type = "real",
  read = list(
    list(
      tag = "CD4Distribution15_49 MV2",
      dims = list("neg_adult_cd4_categories", "g_both_first", "art_status", "years"),
      skip = list(
        rows =  c(25,50)
      )
    ),
    list(
      tag = "CD4Distribution15_49 MV",
      dims = list("cd4_count_dist_old",  "years"),
      skip = list(
        rows =  c(9:11, 19:28, 37:39,
                  47:56, 65:67, 75:84)
      )
    )
  )
),
 aids_deaths_no_art_single_age = list(
   ##Note: not available before 2022
   allow_null = TRUE,
   type = "real",
   read = list(
     list(
       tag = "AIDSDeathsNoARTSingleAge MV",
       dims = list("a_80plus_no_80", "g_both_first", "years"),
       skip = list(
         rows = c(82,163)
       ),
       start_offset = list(row = 2)
     )
   )
 ),
aids_deaths_art_single_age = list(
  ##Note: not available before 2022
  allow_null = TRUE,
  type = "real",
  read = list(
    list(
      tag = "AIDSDeathsARTSingleAge MV",
      dims = list("a_80plus_no_80", "g_both_first", "years"),
      skip = list(
        rows = c(82,163)
      ),
      start_offset = list(row = 2)
    )
  )
),
child_art_calc = list(
  type = "real",
  read = list(
    list(
      tag = "ChildARTCalc MV2",
      dims = list("g_both_first", "child_art_cats", "years"),
      start_offset = list(row = 1)
    ),
    list(
      tag = "ChildARTCalc MV",
      dims = list("g_both_first", "child_art_cats", "years"),
      skip = list(
        rows = 1 + c(dim_vars$g_both_first$length + 1) * 0:3
      ),
      start_offset = list(row = 1)
    )
  )
)
)
}

get_years_cfg <- function() {
  list(
    first_year = list(
      type = "int",
      read = list(
        list(tag = "FirstYear MV", start_offset = list(row = 1)),
        list(tag = "FirstYear MV2")
      )
    ),
    final_year = list(
      type = "int",
      read = list(
        list(tag = "FinalYear MV", start_offset = list(row = 1)),
        list(tag = "FinalYear MV2")
      )
    )
  )
}

get_static_dim_vars <- function() {
  adult_cd4_categories = c("500+", "350-500", "250-349", "200-249", "100-199", "50-99", "50-")
  pmtct_options = c("Single dose nevirapine", "WHO 2006 dual ARV regimen",
    "Option A", "Option B",
    "ART: Started before pregnancy",
    "ART: Started during pregnancy >4 weeks",
    "ART: Started during pregnancy <4 weeks")
  both_sex_vec <- c("both", "male", "female")

  list(
    a = list(length = 81, labels = as.character(0:80)),
    a_80plus = list(length = 82, labels = c(0:80, "81+")),
    a_80plus_no_80 = list(length = 81, labels = c(0:79, "80+")),
    a_15to49_5year = list(length = 7, labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
    a_15to49_custom = list(length = 6, labels = c("15-17", "18-19", "20-24", "25-29", "30-34", "35-49")),
    a_0to80plus_5year = list(length = 17, labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+")),
    a_15to35_custom = list(length = 6, labels = c("15-17", "18-19", "20-24", "25-29", "30-34", "35-49")),
    a_15to45plus_10year = list(length = 4, labels = c("15-24", "25-34", "35-44", "45+")),
    a_0to14 = list(length = 15, labels = 0:14),
    a_0to14_coarse = list(length = 3, labels = c("00-04", "05-09", "10-14")),
    a_0to14_coarse_2 = list(length = 3, labels = c("00-02", "03-04", "05-14")),
    g = list(length = 2, labels = c("male", "female")),
    g_both_middle = list(length = 3, labels = c("male", "both", "female")),
    g_both_first = list(length = 3, labels = both_sex_vec),
    incid_type = list(length = 6, labels = list("Direct incidence input (15 - 49)", "EPP Incidence (15 - 49)", "AEM Incidence (15 - 49)", "Fit incidence to CSAVR", "Fit incidence to mortality data", "ECDC")),
    alloc_type = list(length = 2, labels = list("Expected mortality", "Proportion eligible")),
    cd4_count = list(length = 7, labels = adult_cd4_categories),
    neg_adult_cd4_categories = list(length = 8, labels = c("HIV neg", adult_cd4_categories)),
    cd4_count_dist_old = list(length = 45, labels = paste0(rep(c(rep("HIV, ", 8), rep("ART, ", 7)), 3),
                                                           rep(paste0(both_sex_vec, ": "), each = 15),
                                                       rep(c("HIV neg", rep(adult_cd4_categories, 2)), 3))),
    cd4_prog = list(length = 6, labels = c("500+ -> 350-500", "350-500 -> 250-349", "250-349 -> 200-249", "200-249 -> 100-199", "100-199 -> 50-99", "50-99 -> 50-")),
    cd4_perc_0to4 = list(length = 7, labels = c(">30", "26-30", "21-25", "16-20", "11-15", "5-10", "<5")),
    cd4_prog_0to14 = list(length = 17, labels = paste(c(rep("Age: 0-2, CD4 prog", 6), rep("Age: 3-4, CD4 prog", 6),
                                                       rep("Age: 5-14, CD4 prog", 5)),
                                                      c(rep(c(">30 -> 26-30", "26-30 -> 21-25", "21-25 -> 16-20", "16-20 -> 11-15", "11-15 -> 5-10", "5-10 -> <5"), 2),
                                                           c("1000+ -> 750-999", "750-999 -> 500-749", "500-749 -> 350-499", "350-499 -> 200-349", "200-349 -> <200")))),
    cd4_prog_0to14_old = list(length = 20, labels = paste(c(rep("Age: 0-2, CD4 prog", 7), rep("Age: 3-4, CD4 prog", 7),
                                                        rep("Age: 5-14, CD4 prog", 6)),
                                                      c(rep(c(">30 -> 26-30", "26-30 -> 21-25", "21-25 -> 16-20", "16-20 -> 11-15", "11-15 -> 5-10", "5-10 -> <5", "final"), 2),
                                                        c("1000+ -> 750-999", "750-999 -> 500-749", "500-749 -> 350-499", "350-499 -> 200-349", "200-349 -> <200", "final")))),
    art_dur_2 = list(length = 2, labels = c("0-12mos", "12mos+")),
    pops_for_treat = list(length = 7, labels = c("Pregnant women", "TB/HIV co-infected", "Discordant couples", "Sex workers", "Men who have sex with men", "Injecting drug users", "Other population")),
    elig_perc_year = list(length = 3, labels = c("eligibility", "percent", "year")),
    vt_trt = list(length = 11, labels = c("CD4 <200", "CD4 200-350", "CD4 >350", "Incident infections", "Single dose nevirapine", "WHO 2006 dual ARV regimen", "Option A", "Option B", "ART: Started before pregnancy", "ART: Started during pregnancy >4 weeks", "ART: Started during pregnancy <4 weeks")),
    vt_trt_cat = list(length = 3, labels = c("Perinatal", "Breastfeeding (per month) <350", "Breastfeeding (per month) >=350")),
    vt_time = list(length = 4, labels = c("Perinatal", "Breastfeeding, <6MOS", "Breastfeeding, 7-12MOS", "Breastfeeding, >12MOS")),
    vt_time_old = list(length = 3, labels = c("Perinatal", "Breastfeeding, <12MOS", "Breastfeeding, >12MOS")), ##Note: this is a guess on what this stands for, unable to open that file
    pmtct_editor_order = list(length = 26, labels = c("No prophylaxis- Percent",
                                             paste0(rep(pmtct_options, each = 2), rep(c("- Number", "- Percent"), 7)), "Total prophlyaxis- Number",
                                                      "Postnatal: No prophylaxis- Percent", paste0(rep("Postnatal ", 2), rep(c("Option A", "Option B"), each = 2), rep(c("- Number", "- Percent"), 2)), "Total postnatal prophlyaxis- Number",
                                             paste0("Monthly dropout breastfeeding: ",
                                                    c("Option A", "Option B",
                                                      "ART 0-12 months breastfeeding",
                                                      "ART 12+ months breastfeeding")))),
    ##TODO: Flag with Jeff, assume the same for duration of breastfeeding? modify in the input script?
    pmtct_editor_order_old = list(length = 25, labels = c("No prophylaxis- Percent",
                                                      paste0(pmtct_options, rep(c("- Number", "- Percent"), 7)), "Total prophlyaxis- Number",
                                                      "Postnatal: No prophylaxis- Percent", paste0(c("Option A", "Option B"), rep(c("- Number", "- Percent"), 2)), "Total postnatal prophlyaxis- Number",
                                                      paste0("Monthly dropout breastfeeding: ",
                                                             c("Option A", "Option B",
                                                               "ART")))),
    pmtct_retained_delivery = list(length = 2, labels = c("Percent already on ART retained at delivery",
                                                          "Percent starting ART retained at delivery")),
    bf_months = list(length = 18, labels = paste0("VT_MOS_", c("00_01", "02_03", "04_05",
                                                               "06_07", "08_09", "10_11",
                                                               "12_13", "14_15", "16_17",
                                                               "18_19", "20_21", "22_23",
                                                               "24_25", "26_27", "28_29",
                                                               "30_31", "32_33", "34_35"))),
    paed_treatment = list(length = 5, labels = c("Cotrim", "ART: 0-14y", "ART: 0-4y", "ART: 5-9y", "ART: 10-14y")),
    paed_treatment_old = list(length = 2, labels = c("Cotrim", "ART: 0-14y")),
    paed_art_mort = list(length = 4, labels = c("Age <5, <12 months on ART",
                                                "Age <5, 12+ months on ART",
                                                "Age >=5, <12 months on ART",
                                                "Age >=5, 12+ months on ART")),
   paed_ordinal_cd4_cat = list(length = 7, labels = 1:7),
   paed_ordinal_cd4_cat_old = list(length = 21, labels = rep(1:7,3)),
   paed_age_time_art_cd4 = list(length = 33, labels = paste0(
                                        c(rep(c("0: ", "1-2: ", "3-4: "), each = 7), rep(c("5-9: ", "10-14: "), each = 6)),
                                        c(rep(c(">30", "26-30", "21-25", "16-20", "11-15", "5-10", "<5"),3),
                                          rep(c(">1000", "750-999", "500-749", "350-499", "200-349", "<200"),2)))),
   ##leaving like this as its easier to access in the leapfrog code
   paed_cd4_art_elig = list(length = 8, labels = c(paste0(
      rep(c("CD4 count: ", "CD4 percent: "), each = 4),
      rep(c("<11 MOS", "12-35 MOS", "35-39 MOS", ">=5 YRS"), 2)))),
   cotrim_years_effective = list(length = 5, labels = paste0("Year ", 1:5)),
   art_status = list(length = 2, labels = c("no art", "art")),
   age_by_sex_both = list(length = 243,labels =  paste0(rep(c(0:79, "80+"), each = 3), rep(paste0(" ", both_sex_vec), 81))),
   age_by_sex_both1 = list(length = 243,labels =  paste0(rep(c(0:79, "80+"), 3), rep(paste0(" ", both_sex_vec), each = 81))),
   age_by_sex = list(length = 162,labels =  paste0(rep(c(0:79, "80+"), each = 2), rep(c(" male", " female"), 81))),
   child_art_cats = list(length = 4, labels = c("Children needing cotrim (0-14): ",
                                                            "Children receiving cotrim (0-14): ",
                                                            "Children needing ART (0-14): ",
                                                            "Children receiving ART (0-14): "))
  )
}

process_pjnz <- function(pjnz, use_coarse_age_groups = FALSE) {
  dp <- read_dp(pjnz)
  dat <- parse_dp(dp)
  dim_vars <- dat$dim_vars

  pars <- lapply(dat$data, function(x) if (is.null(x)) NULL else x$data)
  names(pars) <- names(dat$data)

  pars$projection_start_year <- pars$first_year
  # BigPop MV2 has both stratification as well
  pars$basepop <- pars$big_pop[, c("male", "female"), ]
  pars$Sx <- pars$surv_rate

  proportion_of_male_births <- pars$sex_birth_ratio / (pars$sex_birth_ratio + 100)
  pars$births_sex_prop <- rbind(
    male = proportion_of_male_births,
    female = 1 - proportion_of_male_births
  )

  asfr <- pars$asfr / 100
  asfr <- apply(asfr / 5, 2, rep, each = 5)
  asfr_sum <- colSums(asfr)
  asfr_sum[asfr_sum == 0] <- 1
  pars$asfr <- sweep(asfr, 2, pars$tfr / asfr_sum, "*")
  dimnames(pars$asfr) <- list(age = 15:49, year = dim_vars$years$labels)

  pars$migr_age_dist <- pars$migr_age_dist / 100
  migr_age_dist_sum <- colSums(pars$migr_age_dist)
  migr_age_dist_sum[migr_age_dist_sum == 0] <- 1
  pars$migr_age_dist <- sweep(pars$migr_age_dist, 2:3, migr_age_dist_sum, "/")

  netmigr_5year <- sweep(pars$migr_age_dist, 2:3, pars$migr_rate, "*")
  netmigr <- array(
    dim = c(dim_vars$a_80plus_no_80$length, dim_vars$g$length, dim_vars$years$length),
    dimnames = list(
      age = dim_vars$a_80plus_no_80$labels,
      sex = dim_vars$g$labels,
      year = dim_vars$years$labels
    )
  )
  netmigr[1:80, , ] <- apply(netmigr_5year[1:16, , ], 2:3, beers::beers_sub_ordinary)
  netmigr[81, , ] <- netmigr_5year[17, , ]

  u5prop <- array(dim = c(5, 2))
  u5prop[1, ] <- pars$Sx[1, , 1] * 2
  u5prop[2, ] <- pars$Sx[2, , 1] * u5prop[1, ]
  u5prop[3, ] <- pars$Sx[3, , 1] * u5prop[2, ]
  u5prop[4, ] <- pars$Sx[4, , 1] * u5prop[3, ]
  u5prop[5, ] <- pars$Sx[5, , 1] * u5prop[4, ]
  u5prop <- sweep(u5prop, 2, colSums(u5prop), "/")

  netmigr[1:5, 1, ] <- u5prop[, 1, drop = FALSE] %*% netmigr_5year[1, 1, ]
  netmigr[1:5, 2, ] <- u5prop[, 2, drop = FALSE] %*% netmigr_5year[1, 2, ]

  pars$projection_period <- if (pars$version_num > "6.2") "calendar" else "midyear"

  if (pars$projection_period == "midyear") {
    netmigr_adj <- netmigr
    netmigr_adj[-1, , ] <- (netmigr[-1, , ] + netmigr[-81, , ]) / 2
    netmigr_adj[1, , ] <- netmigr[1, , ] / 2
    netmigr_adj[81, , ] <- netmigr_adj[81, , ] + netmigr[81, , ] / 2
    pars$netmigr <- netmigr_adj
  } else {
    pars$netmigr <- netmigr
  }

  pars$netmigr_adj <- pars$netmigr

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

  # TODO !!!!!
  projp <- eppasm::read_hivproj_param(pjnz)
  pars$incrr_sex <- projp$incrr_sex

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
    1, c(3, dim_vars$years$length), list(
      artdur = c("ART0MOS", "ART6MOS", "ART1YR"),
      year = dim_vars$years$labels
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
    pars$cd4_nonaids_excess_mort[, , "male"] <- aperm(
      pars$adult_non_aids_excess_mort[1, , , ], c(2, 3, 1)
    )
    pars$cd4_nonaids_excess_mort[, , "female"] <- aperm(
      pars$adult_non_aids_excess_mort[3, , , ], c(2, 3, 1)
    )
    pars$art_nonaids_excess_mort[, , "male"] <- aperm(
      pars$adult_non_aids_excess_mort[2, , , ], c(2, 3, 1)
    )
    pars$art_nonaids_excess_mort[, , "female"] <- aperm(
      pars$adult_non_aids_excess_mort[4, , , ], c(2, 3, 1)
    )
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

  ## Add in pediatric components
  frr_agecat <- seq(15, 45, 5)
  h_fert_idx <- which((15L - 1 + cumsum(h_ag_span_coarse)) %in% 15:49)
  age_band_width <- length(15:49) / length(frr_agecat)
  fert_rat_h_ag <- findInterval(
    15L + cumsum(h_ag_span_coarse[h_fert_idx]) - h_ag_span_coarse[h_fert_idx],
    frr_agecat
  )
  pars$fert_rat_coarse <- array(1, c(length(h_fert_idx), dim_vars$years$length))
  pars$fert_rat_coarse[, ] <- pars$hivtfr[fert_rat_h_ag, ]
  pars$fert_rat_full <- apply(pars$hivtfr, 2, rep, each = age_band_width)

  if (is.null(pars$fert_cd4_discount)) {
    pars$cd4fert_rat <- rep(1, dim_vars$cd4_count$length)
  } else {
    pars$cd4fert_rat <- pars$fert_cd4_discount
  }

  if (is.null(pars$ratio_women_on_art)) {
    pars$ratio_women_on_art <- rep(1, dim_vars$a_15to49_5year$length)
  } else if (dat$data$ratio_women_on_art$tag == "RatioWomenOnART MV") {
    pars$ratio_women_on_art <- rep(pars$ratio_women_on_art, dim_vars$a_15to49_5year$length)
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
