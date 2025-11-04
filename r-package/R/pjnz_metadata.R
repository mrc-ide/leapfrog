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
  metadata <- c(get_years_cfg(), get_pars_metadata(dim_vars))

  ret <- lapply(names(metadata), function(name) {
    get_data_from_cfg(name, metadata[[name]], dim_vars, dp)
  })
  names(ret) <- names(metadata)
  list(data = ret, dim_vars = dim_vars)
}

get_pars_metadata <- function(dim_vars) {
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
          column_dims = list("cd4_prog", "a_15to45plus_10year")
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
          dims = list("art_off_on", "g", "cd4_count", "a_15to45plus_10year"),
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
  list(
    a = list(length = 81, labels = as.character(0:80)),
    a_80plus = list(length = 82, labels = c(0:80, "81+")),
    a_80plus_no_80 = list(length = 81, labels = c(0:79, "80+")),
    a_15to49_5year = list(length = 7, labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
    a_15to49_custom = list(length = 6, labels = c("15-17", "18-19", "20-24", "25-29", "30-34", "35-49")),
    a_0to80plus_5year = list(length = 17, labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+")),
    a_15to35_custom = list(length = 6, labels = c("15-17", "18-19", "20-24", "25-29", "30-34", "35-49")),
    a_15to45plus_10year = list(length = 4, labels = c("15-24", "25-34", "35-44", "45+")),
    g = list(length = 2, labels = c("male", "female")),
    g_both_middle = list(length = 3, labels = c("male", "both", "female")),
    g_both_first = list(length = 3, labels = c("both", "male", "female")),
    cd4_count = list(length = 7, labels = c("500+", "350-500", "250-349", "200-249", "100-199", "50-99", "50-")),
    cd4_prog = list(length = 6, labels = c("500+ -> 350-500", "350-500 -> 250-349", "250-349 -> 200-249", "200-249 -> 100-199", "100-199 -> 50-99", "50-99 -> 50-")),
    art_dur_2 = list(length = 2, labels = c("0-12mos", "12mos+")),
    art_off_on = list(length = 2, labels = c("off art", "on art")),
    pops_for_treat = list(length = 7, labels = c("Pregnant women", "TB/HIV co-infected", "Discordant couples", "Sex workers", "Men who have sex with men", "Injecting drug users", "Other population")),
    elig_perc_year = list(length = 3, labels = c("eligibility", "percent", "year"))
  )
}

process_pjnz <- function(pjnz) {
  dp <- read_dp(pjnz)
  dat <- parse_dp(dp)
  dim_vars <- dat$dim_vars

  pars <- lapply(dat$data, function(x) if (is.null(x)) NULL else x$data)
  names(pars) <- names(dat$data)

  pars$projection_start_year <- pars$first_year
  # BigPop MV2 has both stratification as well
  pars$base_pop <- pars$big_pop[, c("male", "female"), ]
  pars$Sx <- pars$surv_rate

  proporiion_of_male_births <- pars$sex_birth_ratio / (pars$sex_birth_ratio + 100)
  pars$births_sex_prop <- rbind(
    male = proporiion_of_male_births,
    female = 1 - proporiion_of_male_births
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
  pars$test_nm5 <- netmigr_5year
  netmigr[1:80, , ] <- apply(netmigr_5year[1:16, , ], 2:3, beers::beers_sub_ordinary)
  netmigr[81, , ] <- netmigr_5year[17, , ]
  pars$test_nm <- netmigr

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

  if (demp$projection_period == "midyear") {
    netmigr_adj <- netmigr
    netmigr_adj[-1, , ] <- (netmigr[-1, , ] + netmigr[-81, , ]) / 2
    netmigr_adj[1, , ] <- netmigr[1, , ] / 2
    netmigr_adj[81, , ] <- netmigr_adj[81, , ] + netmigr[81, , ] / 2
    pars$netmigr <- netmigr_adj
  } else {
    pars$netmigr <- netmigr
  }

  pars$netmigr_adj <- pars$netmigr

  pars
}
