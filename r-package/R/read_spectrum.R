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

#' Get DP data
#'
#' This function operates as a pass-through if the file is already
#' a parsed DP file, or extracts it from the CSV if it is a filepath
#' to a PJNZ file.
#'
#' @param dp_pjnz either a `"spectrum_dp"` object or a file path to
#'   a PJNZ file.
#'
#' @return An object of class `"spectrum_dp"`, which is the .DP file
#'   read by `read.csv()`.
#'
#'
#' @noRd
get_dp_data <- function(dp_pjnz) {

  if (inherits(dp_pjnz, "spectrum_dp")) {
    dp_pjnz
  } else if (is.character(dp_pjnz) && file.exists(dp_pjnz)) {
    read_dp(dp_pjnz)
  } else {
    stop("Valid input not found. Provide path to PJNZ file.")
  }

}

exists_dptag <- function(dp, tag, tagcol = 1) {

  stopifnot(inherits(dp, "spectrum_dp"))
  stopifnot(is.character(tag))

  tag %in% dp[, tagcol]
}

exists_dpdescription <- function(dp, description, descriptioncol = 2) {

  stopifnot(inherits(dp, "spectrum_dp"))
  stopifnot(is.character(description))

  description %in% dp[, descriptioncol]
}

dpsub <- function(dp, tag, rows, cols, tagcol = 1) {

  stopifnot(inherits(dp, "spectrum_dp"))
  stopifnot(is.character(tag))
  all.equal(rows, as.integer(rows))
  all.equal(cols, as.integer(cols))

  dp[which(dp[, tagcol] == tag) + rows, cols]
}

dpdescription <- function(dp, description, rows, cols, tagcol = 2) {

  stopifnot(inherits(dp, "spectrum_dp"))
  all.equal(rows, as.integer(rows))
  all.equal(cols, as.integer(cols))

  dp[which(dp[, tagcol] == description) + rows, cols]
}

get_dp_years <- function(dp.x) {
  yr_start <- as.integer(dpsub(dp = dp.x, tag =  "<FirstYear MV2>",rows = 2,cols = 4))
  yr_end <- as.integer(dpsub(dp = dp.x,tag =  "<FinalYear MV2>", rows = 2,cols = 4))

  proj_years <- yr_start:yr_end
  time_data_idx <- 4 + 1:length(proj_years) - 1

  list(proj_years = proj_years,
       time_data_idx = time_data_idx)
}


#' Read Spectrum programme data inputs
#'
#' @param dp An `"spectrum_dp"` object created by `read_dp()` or a path to a Spectrum
#'   PJNZ file.
#'
#' @return
#'
#' Most functions return a matrix in which rows are each indicator and columns are years
#'
#' `dp_read_abortion()` returns a list in which the vector `pregtermabortion` is the
#' numerical values and `pregtermabortion_isperc` is a logical indicating whether each
#' value corresponds to a number (count) input or percentage input.
#'
#' @details
#'
#' These functions accepts either a path to a PJNZ file or an already parsed `.DP`
#' file read using [`read_dp()`] with class `"spectrum_dp"`. Unzipping and reading the
#' `.DP` file is the slow part of the function, so if reading lots of data from
#' the same file, it will be most efficient to read it once and pass that to the
#' functions.
#'
#' @examples \dontrun{
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "frogger")
#' dp <- read_dp(pjnz)
#' dp_anc_testing <- dp_read_anc_testing(dp)
#' dp_pmtct <- dp_read_pmtct(dp)
#' dp_read_pmtct <- dp_read_pmtct_retained(dp)
#' dp_abortion <- dp_read_abortion(dp)
#' dp_notbreastfeeding <- dp_read_breastfeeding(dp)
#' dp_childart <- dp_read_childart(dp)
#' dp_childltfu <- dp_read_childltfu(dp)
#'
#' ## Can either pass PJNZ path or parsed "spectrum_dp" object
#'
#' anc_testing1 <- dp_read_anc_testing(pjnz)
#' anc_testing2 <- dp_read_anc_testing(dp)
#' all.equal(anc_testing1, anc_testing2)
#' }
#' @keywords internal
dp_read_anc_testing <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  anc_indicators <- c("anc_clients", "anc_tested", "anc_tested_pos", "anc_known_pos",
                      "anc_prevalence", "anc_retested", "anc_retested_pos",
                      "facility_births", "anc_known_neg")

  if (exists_dptag(dp, "<ANCTestingValues MV>")) {
    anc_testing <- dpsub(dp, "<ANCTestingValues MV>", 2:5, dpy$time_data_idx)
    anc_testing <- df_as_numeric(anc_testing)
    dimnames(anc_testing) <- list(indicator = anc_indicators[1:4], year = dpy$proj_years)
  } else if (exists_dptag(dp, "<ANCTestingValues MV2>")) {
    anc_testing <- dpsub(dp, "<ANCTestingValues MV2>", 2:5, dpy$time_data_idx)
    anc_testing <- df_as_numeric(anc_testing)
    dimnames(anc_testing) <- list(indicator = anc_indicators[1:4], year = dpy$proj_years)
  } else if (exists_dptag(dp, "<ANCTestingValues MV4>")) {
    anc_testing <- dpsub(dp, "<ANCTestingValues MV4>", 2:10, dpy$time_data_idx)
    anc_testing <- df_as_numeric(anc_testing)
    dimnames(anc_testing) <- list(indicator = anc_indicators, year = dpy$proj_years)
  } else {
    stop("ANC testing tag not recognized. Function probably needs update for this .DP file.")
  }

  anc_testing[anc_testing == -9999] <- NA_real_

  anc_testing
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_pmtct <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  pmtct_indicators <- c("pmtct_noprophylaxis_percent",
                        "pmtct_singledosenvp_number",
                        "pmtct_singledosenvp_percent",
                        "pmtct_dualarv_number",
                        "pmtct_dualarv_percent",
                        "pmtct_optiona_number",
                        "pmtct_optiona_percent",
                        "pmtct_optionb_number",
                        "pmtct_optionb_percent",
                        "pmtct_alreadyart_number",
                        "pmtct_alreadyart_percent",
                        "pmtct_artbefore4weeks_number",
                        "pmtct_artbefore4weeks_percent",
                        "pmtct_artafter4weeks_number",
                        "pmtct_artafter4weeks_percent",
                        "pmtct_total_number",
                        "pmtct_postnatal_noprophylaxis_percent",
                        "pmtct_postnatal_optiona_number",
                        "pmtct_postnatal_optiona_percent",
                        "pmtct_postnatal_optionb_number",
                        "pmtct_postnatal_optionb_percent",
                        "pmtct_postnatal_total_number",
                        "pmtct_postnatal_monthlydropout_optiona",
                        "pmtct_postnatal_monthlydropout_optionb",
                        "pmtct_postnatal_monthlydropout_art0to12months",
                        "pmtct_postnatal_monthlydropout_art12plusmonths")


  ## Note: these values start 1 column later than other arrays in the .DP file
  ## If value is 0, interpret as not entered (NA)

  if (exists_dptag(dp, "<ARVRegimen MV3>")) {
    pmtct_arv <- dpsub(dp, "<ARVRegimen MV3>", 2:27, dpy$time_data_idx + 1L)
    pmtct_arv <- df_as_numeric(pmtct_arv)
    dimnames(pmtct_arv) <- list(indicator = pmtct_indicators, year = dpy$proj_years)
  } else {
    stop("PMTCT ARV tag not recognized. Function probably needs update for this .DP file.")
  }

  pmtct_arv[pmtct_arv == 0.0] <- NA_real_

  pmtct_arv
}


#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_pmtct_retained <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  indicator_names <- c("pmtct_retained_alreadyart", "pmtct_retained_newart")

  if (exists_dptag(dp, "<PercentARTDelivery MV>")) {
    pmtct_retained <- dpsub(dp, "<PercentARTDelivery MV>", 2:3, dpy$time_data_idx)
  } else {
    stop("PMTCT retained at delivery tag not recognized. Function probably needs update for this .DP file.")
  }

  pmtct_retained <- df_as_numeric(pmtct_retained)
  dimnames(pmtct_retained) <- list(indicator = indicator_names, year = dpy$proj_years)

  pmtct_retained
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_abortion <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  ## Note: If 0s are entered for both, output = 1 defaults to "percentage"
  if (exists_dptag(dp, "<PregTermAbortionPerNum MV2>")) {
    pregtermabortion_ispercent <- dpsub(dp, "<PregTermAbortionPerNum MV2>", 2, dpy$time_data_idx)
  } else {
    stop("PregTermAbortionPerNum tag not found. Function probably needs update for this .DP file.")
  }

  if (exists_dptag(dp, "<PregTermAbortion MV3>")) {
    pregtermabortion <- dpsub(dp, "<PregTermAbortion MV3>", 2, dpy$time_data_idx)
  } else {
    stop("PregTermAbortionPerNum tag not found. Function probably needs update for this .DP file.")
  }

  pregtermabortion <- stats::setNames(as.numeric(pregtermabortion), dpy$proj_years)
  pregtermabortion_ispercent <- stats::setNames(as.logical(pregtermabortion_ispercent), dpy$proj_years)

  list(pregtermabortion = pregtermabortion,
       pregtermabortion_ispercent = pregtermabortion_ispercent)
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_mothers_reallocated <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)


  if (exists_dptag(dp, "<DP_TGX_PatientsReallocated_MV>")) {
    patients_reallocated <- dpsub(dp, "<DP_TGX_PatientsReallocated_MV>", 2, dpy$time_data_idx)
  } else {
    stop("PatientsReallocated tag not found. Function probably needs update for this .DP file.")
  }

  patients_reallocated <- stats::setNames(as.numeric(patients_reallocated), dpy$proj_years)

  patients_reallocated
}


#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_breastfeeding <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  agemonth_cat <- c("M00_01", "M02_03", "M04_05", "M06_07", "M08_09", "M10_11",
                    "M12_13", "M14_15", "M16_17", "M18_19", "M20_21", "M22_23",
                    "M24_25", "M26_27", "M28_29", "M30_31", "M32_33", "M34_35")

  ## Note: If 0s are entered for both, output = 1 defaults to "percentage"
  if (exists_dptag(dp, "<InfantFeedingOptions MV>")) {
    notbreastfeeding_percent <- dpsub(dp, "<InfantFeedingOptions MV>", 2:37, dpy$time_data_idx)

  } else {
    stop("Not Breastfeeding Percent tag not found. Function probably needs update for this .DP file.")
  }

  notbreastfeeding_percent <- df_as_numeric(notbreastfeeding_percent)
  notbreastfeeding_percent_noarv <- notbreastfeeding_percent[1:18, ]
  notbreastfeeding_percent_arv <- notbreastfeeding_percent[18 + 1:18, ]

  dn <- list(child_age_months = agemonth_cat, year = dpy$proj_years)
  dimnames(notbreastfeeding_percent_noarv) <- dn
  dimnames(notbreastfeeding_percent_arv) <- dn

  list(notbreastfeeding_percent_noarv = notbreastfeeding_percent_noarv,
       notbreastfeeding_percent_arv = notbreastfeeding_percent_arv)
}


#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_childart <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  indicator_names <- c("childart_ctx",
                       "childart_art0to14",
                       "childart_art0to4",
                       "childart_art5to9",
                       "childart_art10to14")

  if (exists_dptag(dp, "<ChildTreatInputs MV3>")) {
    childart <- dpsub(dp, "<ChildTreatInputs MV3>", 2:6, dpy$time_data_idx)
  } else {
    stop("Child ART input tag not recognized. Function probably needs update for this .DP file.")
  }

  childart <- df_as_numeric(childart)
  dimnames(childart) <- list(indicator = indicator_names, year = dpy$proj_years)
  childart[childart == -9999] <- NA_real_

  if (exists_dptag(dp, "<ChildARTByAgeGroupPerNum MV2>")) {
    childart_ispercent <- dpsub(dp, "<ChildARTByAgeGroupPerNum MV2>", 2:6, dpy$time_data_idx)
  } else {
    stop("Child ART input tag not recognized. Function probably needs update for this .DP file.")
  }

  childart_ispercent <- df_as_logical(childart_ispercent)
  dimnames(childart_ispercent) <- list(indicator = indicator_names, year = dpy$proj_years)

  list(childart = childart,
       childart_ispercent = childart_ispercent)
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_childltfu <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  if (exists_dptag(dp, "<PercLostFollowupChild MV>")) {
    childart_ltfu <- dpsub(dp, "<PercLostFollowupChild MV>", 2, dpy$time_data_idx)
  } else {
    stop("Child ART LTFU tag not recognized. Function probably needs update for this .DP file.")
  }


  childart_ltfu <- as.numeric(childart_ltfu)
  names(childart_ltfu) <- dpy$proj_years

  childart_ltfu
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_art_dist <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  if (exists_dptag(dp, "<ChildARTDist MV>")) {
    child_mort_mult <- dpsub(dp, "<ChildARTDist MV>", 2:16, dpy$time_data_idx)
    child_mort_mult <- df_as_numeric(child_mort_mult)
    dimnames(child_mort_mult) <- list(age = 0:14, year = dpy$proj_years)
  } else {
    stop("Child ART distribution tag not recognized. Function probably needs update for this .DP file.")
  }

  child_mort_mult
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_child_mort_mult <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  indicator_names <- c('Age <5, <12 months on ART',
                       'Age <5, 12+ months on ART',
                       'Age >=5, <12 months on ART',
                       'Age >=5, 12+ months on ART')

  if (exists_dptag(dp, "<ChildMortalityRates MV2>")) {
    child_mort_mult <- dpsub(dp, "<ChildMortalityRates MV2>", 2:5, dpy$time_data_idx)
    child_mort_mult <- df_as_numeric(child_mort_mult)
    dimnames(child_mort_mult) <- list(mrr = indicator_names, year = dpy$proj_years)
  } else {
    stop("Child ART multiplier tag not recognized. Function probably needs update for this .DP file.")
  }

  child_mort_mult
}


#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_nosocom_infections <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)


  if (exists_dptag(dp, "<NosocomialInfectionsByAge MV>")) {
    ##only extracting 0-4 for right now
    nosocomial_inf <- dpsub(dp, "<NosocomialInfectionsByAge MV>", 2, dpy$time_data_idx)
    nosocomial_inf <- df_as_numeric(nosocomial_inf)
    names(nosocomial_inf) <- dpy$proj_years
  } else {
    stop("Nosocomial infections tag not recognized. Function probably needs update for this .DP file.")
  }

  nosocomial_inf
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_mtct_rates <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)


  if (exists_dpdescription(dp, "Peripartum and breastfeeding transmission rates (%)")) {
    ##only extracting 0-4 for right now
    mtct_rates <- dpdescription(dp, "Peripartum and breastfeeding transmission rates (%)" , 1:11, 4:6)
    mtct_rates <- df_as_numeric(mtct_rates)
  } else {
    stop("MTCT rates description not recognized. Function probably needs update for this .DP file.")
  }

  mtct_rates
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_paed_cd4_dist <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)


  if (exists_dpdescription(dp, "Distribution of new infections by CD4 percent for Children")) {
    ##only extracting 0-4 for right now
    cd4_dist <- dpdescription(dp, "Distribution of new infections by CD4 percent for Children" , 1, 4:10)
    cd4_dist <- df_as_numeric(cd4_dist)
  } else if(exists_dpdescription(dp, "Répartition des nouvelles infections par pourcentage de CD4 chez les enfants")) {
    cd4_dist <- dpdescription(dp, "Répartition des nouvelles infections par pourcentage de CD4 chez les enfants", 1, 4:10)
    cd4_dist <- df_as_numeric(cd4_dist)
  } else {
    stop("CD4 distribution for paeds description not recognized. Function probably needs update for this .DP file.")
  }

  cd4_dist
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_paed_cd4_prog <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  cd4 <- c(rep(c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10'), 2), NA,
           c('>1000', '750-999', '500-749', '350-499', '200-349'))
  ages <- c('<5', '5-14')

  if (exists_dpdescription(dp, "Annual rate of progression to next lower CD4 category for Children")) {
    ##only extracting 0-4 for right now
    cd4_prog <- dpdescription(dp, "Annual rate of progression to next lower CD4 category for Children", 2:3, 4:21)
    cd4_prog <- df_as_numeric(cd4_prog)
    dimnames(cd4_prog) <- list(sex = c('Male', 'Female'), cd4_cat = cd4)


  } else {
    stop("CD4 distribution for paeds description not recognized. Function probably needs update for this .DP file.")
  }

  cd4_prog
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_paed_cd4_mort <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  cd4 <- 1:7
  trans_type = rep(c('perinatal', 'bf 0-6', 'bf 7-12', 'bf 12+'),3)

  if (exists_dpdescription(dp, "Annual probability of HIV-related mortality among those not on ART by CD4 category for Children")) {
    ##only extracting 0-4 for right now
    cd4_mort <- dpdescription(dp, "Annual probability of HIV-related mortality among those not on ART by CD4 category for Children", 1:12, 4:10)
    cd4_mort <- df_as_numeric(cd4_mort)
    dimnames(cd4_mort) <- list(trans = trans_type, cd4_cat = cd4)


  } else {
    stop("CD4 mortality for paeds description not recognized. Function probably needs update for this .DP file.")
  }

  cd4_mort
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_paed_art_mort <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  cd4 <- 1:7
  trans_type = rep(c('perinatal', 'bf 0-6', 'bf 7-12', 'bf 12+'),3)

  if (exists_dptag(dp, "<ChildMortByCD4WithART0to6 MV2>")) {
    ##only extracting 0-4 for right now
    art_mort <- dpsub(dp,"<ChildMortByCD4WithART0to6 MV2>",
                              2:3, 4:25)
    art_mort <- df_as_numeric(art_mort)
    males_lt6mo <- array(art_mort[1,], dim = c(7, 3), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4')))
    females_lt6mo <- array(art_mort[2,], dim = c(7, 3), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4')))
    art_mort_lt6mo <- array(0, dim = c(7, 3, 2), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4'), sex = c('Male', 'Female')))
    art_mort_lt6mo[,,1] <- males_lt6mo
    art_mort_lt6mo[,,2] <- females_lt6mo


    ##only extracting 0-4 for right now
    art_mort <- dpsub(dp,"<ChildMortByCD4WithART0to6 MV2>",
                      2:3, 26:38)
    art_mort <- df_as_numeric(art_mort)
    art_mort <- art_mort[,-7]
    adol_males_lt6mo <- array(art_mort[1,], dim = c(6, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14')))
    adol_females_lt6mo <- array(art_mort[2,], dim = c(6, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14')))
    adol_art_mort_lt6mo <- array(0, dim = c(6, 2, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14'), sex = c('Male', 'Female')))
    adol_art_mort_lt6mo[,,1] <- adol_males_lt6mo
    adol_art_mort_lt6mo[,,2] <- adol_females_lt6mo

  } else {
    stop("ART mortality less than 6 months for paeds description not recognized. Function probably needs update for this .DP file.")
  }

  if (exists_dptag(dp, "<ChildMortByCD4WithART7to12 MV>")) {
    ##only extracting 0-4 for right now
    art_mort <- dpsub(dp,"<ChildMortByCD4WithART7to12 MV>",
                      3:4, 4:25)
    art_mort <- df_as_numeric(art_mort)
    males_6to12mo <- array(art_mort[1,], dim = c(7, 3), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4')))
    females_6to12mo <- array(art_mort[2,], dim = c(7, 3), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4')))
    art_mort_6to12mo <- array(0, dim = c(7, 3, 2), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4'), sex = c('Male', 'Female')))
    art_mort_6to12mo[,,1] <- males_6to12mo
    art_mort_6to12mo[,,2] <- females_6to12mo

    art_mort <- dpsub(dp,"<ChildMortByCD4WithART7to12 MV>",
                      3:4, 26:38)
    art_mort <- df_as_numeric(art_mort)
    art_mort <- art_mort[,-7]
    adol_males_6to12mo <- array(art_mort[1,], dim = c(6, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14')))
    adol_females_6to12mo <- array(art_mort[2,], dim = c(6, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14')))
    adol_art_mort_6to12mo <- array(0, dim = c(6, 2, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14'), sex = c('Male', 'Female')))
    adol_art_mort_6to12mo[,,1] <- adol_males_6to12mo
    adol_art_mort_6to12mo[,,2] <- adol_females_6to12mo

  } else {
    stop("ART mortality 6 to 12 months for paeds description not recognized. Function probably needs update for this .DP file.")
  }

  if (exists_dptag(dp, "<ChildMortByCD4WithARTGT12 MV>")) {
    ##only extracting 0-4 for right now
    art_mort <- dpsub(dp,"<ChildMortByCD4WithARTGT12 MV>",
                      3:4, 4:25)
    art_mort <- df_as_numeric(art_mort)
    males_gte12mo <- array(art_mort[1,], dim = c(7, 3), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4')))
    females_gte12mo <- array(art_mort[2,], dim = c(7, 3), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4')))
    art_mort_gte12mo <- array(0, dim = c(7, 3, 2), dimnames = list(cd4 = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'), ages = c('0', '1-2', '3-4'), sex = c('Male', 'Female')))
    art_mort_gte12mo[,,1] <- males_gte12mo
    art_mort_gte12mo[,,2] <- females_gte12mo

    art_mort <- dpsub(dp,"<ChildMortByCD4WithARTGT12 MV>",
                      3:4, 26:38)
    art_mort <- df_as_numeric(art_mort)
    art_mort <- art_mort[,-7]
    adol_males_gte12mo <- array(art_mort[1,],  dim = c(6, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14')))
    adol_females_gte12mo <- array(art_mort[2,], dim = c(6, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14')))
    adol_art_mort_gte12mo <- array(0, dim = c(6, 2, 2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', '<200'), ages = c('5-9', '10-14'), sex = c('Male', 'Female')))
    adol_art_mort_gte12mo[,,1] <- adol_males_gte12mo
    adol_art_mort_gte12mo[,,2] <- adol_females_gte12mo

  } else {
    stop("ART mortality greater than 12 months for paeds description not recognized. Function probably needs update for this .DP file.")
  }

   return(list(hc1_lt6 = art_mort_lt6mo, hc1_6to12 = art_mort_6to12mo, hc1_gte12 = art_mort_gte12mo,
               hc2_lt6 = adol_art_mort_lt6mo, hc2_6to12 = adol_art_mort_6to12mo, hc2_gte12 = adol_art_mort_gte12mo))
}

#' @rdname dp_read_anc_testing
#' @keywords internal
dp_read_paed_art_eligibility <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  specs <- c(paste0('CD4 count: ', c('Age < 11 months', 'Age 12-35 months', 'Age 35-39 months', 'Age >= 5 years')),
             paste0('CD4 percent: ', c('Age < 11 months', 'Age 12-35 months', 'Age 35-39 months', 'Age >= 5 years')))

  if (exists_dpdescription(dp, "Eligibility for treatment - Children")) {
    ##only extracting 0-4 for right now
    art_elig <- dpdescription(dp, "Eligibility for treatment - Children", 1:8, dpy$time_data_idx)
    art_elig <- df_as_numeric(art_elig)
    dimnames(art_elig) <- list(age = specs, years = dpy$proj_years)

    art_elig_age <- dpdescription(dp, "Age below which all HIV+ children should be on treatment (months)", 1, dpy$time_data_idx)
    art_elig_age <- df_as_numeric(art_elig_age)
    names(art_elig_age) <-  dpy$proj_years

  } else {
    stop("CD4 mortality for paeds description not recognized. Function probably needs update for this .DP file.")
  }

  list(cd4_elig = art_elig, age_elig = art_elig_age)
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
#'   demp <- prepare_leapfrog_demp(pjnz)
#'   proj <- prepare_leapfrog_projp(pjnz)
#'   parameters <- c(proj, demp)
#' projp <- prepare_hc_leapfrog_projp(pjnz, parameters)
#' @export
prepare_hc_leapfrog_projp <- function(pjnz, params) {
  dp.x <- get_dp_data(pjnz)

  ## projection parameters
  yr_start <- as.integer(dpsub(dp = dp.x, "<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub(dp = dp.x, "<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  year.idx <- 1:length(proj.years)

  ## Hard coded to expand age groups 15-24, 25-34, 35-44, 45+ to
  ## single-year ages 15:80.
  ## Requires extension for coarse HIV age group stratification
  idx_expand_full <- rep(1:4, times = c(10, 10, 10, 36))
  idx_expand_coarse <- rep(1:4, times = c(3, 2, 2, 2))

  v <- params
  ## paed input
  v$paed_incid_input <- dp_read_nosocom_infections(dp.x)
  v$paed_cd4_dist <- dp_read_paed_cd4_dist(dp.x) / 100

  prog = dp_read_paed_cd4_prog(dp.x)
  paed_cd4_prog <- array(0, dim = c(7,2,2), dimnames = list(cd4pct = c('30plus', '26-30', '21-25', '16-20', '11-15', '5-10', '<5'),
                                                            age = c('0-2', '3-4'),
                                                            sex = c('male', 'female')))
  paed_cd4_prog[,1,1] <- c(prog[1,1:6], 0)
  paed_cd4_prog[,1,2] <- c(prog[2,1:6], 0)

  paed_cd4_prog[,2,1] <- c(prog[1,7:12], 0)
  paed_cd4_prog[,2,2] <- c(prog[2,7:12], 0)

  adol_cd4_prog <- array(c(prog[1,14:18],0), dim = c(6,1,2), dimnames = list(cd4 = c('>1000', '750-999', '500-749', '350-499', '200-349', 'lte200'),
                                                                             age = c('5-14'),
                                                                             sex = c('male', 'female')))

  v$paed_cd4_prog <- paed_cd4_prog
  v$adol_cd4_prog <- adol_cd4_prog

  mort <- dp_read_paed_cd4_mort(dp.x)
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


  mort <- dp_read_paed_art_mort(dp.x)
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



  art_dist_paed <- dp_read_art_dist(dp.x)
  v$art_dist_paed <- art_dist_paed

  ## pull in cotrim coverage numbers
  ctx_pct <- input_childart(dp.x)$ctx_percent
  ctx_pct[is.na(ctx_pct)] <- FALSE
  v$ctx_val_ispercent <- ctx_pct
  v$ctx_val <- input_childart(dp.x)$ctx
  if(any(v$ctx_val_ispercent)){
    v$ctx_val[v$ctx_val_ispercent] <- v$ctx_val[v$ctx_val_ispercent] / 100
  }
  ##cotrim is effective for five years for children not on ART and for four years for children on ART
  ctx_effect <- dpsub(dp = dp.x, "<EffectTreatChild MV>",3:4,4:13)
  off_art_ctx <- sum(as.numeric(unlist(ctx_effect[1,]))) / 5
  on_art_ctx.lte12mo <- sum(as.numeric(unlist(ctx_effect[2,1])))
  on_art_ctx.gte12mo <- sum(as.numeric(unlist(ctx_effect[2,2:5]))) / 4
  ctx_effect <- array(data = c(off_art_ctx, on_art_ctx.lte12mo, on_art_ctx.gte12mo),
                      dim = c(3),
                      dimnames = list(ctx_effect = c('Off ART', 'On ART, lte12mo', 'On ART, gte12mo')))
  v$ctx_effect <- ctx_effect

  ## pull in ART coverage numbers
  art = input_childart(dp.x)
  v$artpaeds_isperc <- art$art_ispercent[2,]
  v$artpaeds_isperc[] <- as.integer(ifelse(art$art_ispercent[2,] == FALSE, 0, 1))
  art$child_art[1,which(art$art_ispercent[2,])] <- art$child_art[1,which(art$art_ispercent[2,])] / 100
  v$paed_art_val <- art$child_art
  v$paed_art_age_spec <- art$age_spec
  v$hc_art_start <- as.integer(unname(which(colSums(art$child_art) > 0)[1]) - 1)

  ##PMTCT
  pmtct_list <- input_pmtct(dp.x)
  pmtct_list <- pmtct_list[c(3,4,1,2,5:7),,]

  if(sum(pmtct_list[,,1]) == 0){
    v$pmtct_input_isperc = rep(F, length(proj.years))
  }else{
    v$pmtct_input_isperc = rep(T, length(proj.years))
  }

  pmtct_new <- array(0, dim = c(7, 61), dimnames = list(pmtct = c("Option A", "Option B", "SDNVP", "Dual ARV", "Option B+: before pregnancy", "Option B+: >4 weeks", "Option B+: <4 weeks")))
  ## pick out which ones were inserted as numbers
  pmtct_new[, which(colSums(pmtct_list)[, 1] > 0)] <- pmtct_list[, (which(colSums(pmtct_list)[, 1] > 0)), 1]
  ## pick out which ones were inserted as percent
  pmtct_new[, which(colSums(pmtct_list)[, 1] == 0)] <- pmtct_list[, which(colSums(pmtct_list)[, 1] == 0), 2]
  v$pmtct <- pmtct_new

  v$pmtct_input_isperc <- !(apply(input_pmtct_ispercent(dp.x), 2, any))

  ##PMTCT dropout
  v$pmtct_dropout <- input_pmtct_retained(dp.x)

  ##rates of MTCT
  mtct_rates_input <- dp_read_mtct_rates(dp.x)
  mtct_trt <- array(data = 0, dim = c(7,7,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'),
                                                              pmtct_reg = c('option A', 'option B', 'single dose nevirapine', 'WHO 2006 dual ARV regimen', 'ART before pregnancy',
                                                                            'ART >4 weeks before delivery', 'ART <4 weeks before delivery'),
                                                              transmission_type = c('perinatal', 'breastfeeding')))
  mtct_trt[,1,1] <-  mtct_rates_input[7,1] /100
  mtct_trt[,2,1] <-  mtct_rates_input[8,1]/100
  mtct_trt[,3,1] <- mtct_rates_input[5,1]/100
  mtct_trt[,4,1] <- mtct_rates_input[6,1]/100
  mtct_trt[3:7,1,2] <- mtct_rates_input[7,3]/100
  mtct_trt[3:7,2,2] <- mtct_rates_input[8,3]/100
  mtct_trt[1:2,3,2] <- mtct_rates_input[5,2]/100
  mtct_trt[3:7,3,2] <- mtct_rates_input[5,3]/100
  mtct_trt[,4,2] <- mtct_rates_input[6,3]/100
  mtct_trt[,5,1] <- mtct_rates_input[9,1]/100
  mtct_trt[,6,1] <- mtct_rates_input[10,1]/100
  mtct_trt[,7,1] <- mtct_rates_input[11,1]/100
  mtct_trt[3:7,5,2] <- mtct_rates_input[9,2]/100
  mtct_trt[3:7,6,2] <- mtct_rates_input[10,2]/100
  mtct_trt[3:7,7,2] <- mtct_rates_input[11,2]/100
  v$pmtct_mtct <- mtct_trt

  mtct <- array(data = NA, dim = c(8,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50', 'INFECTION'), trans_type = c('perinatal', 'bf')))
  mtct[1:7,1] <- c(mtct_rates_input[3,1], mtct_rates_input[3,1],
                   mtct_rates_input[2,1], mtct_rates_input[2,1],
                   rep(mtct_rates_input[1,1],3)) /100
  mtct[1:7,2] <- c(rep(mtct_rates_input[3,3],2), rep(mtct_rates_input[2,2],2), rep(mtct_rates_input[1,2],3)) / 100
  mtct[8,] <- c(mtct_rates_input[4,1], mtct_rates_input[4,2]) / 100
  v$mtct <- mtct

  mort_rr_art <- dp_read_child_mort_mult(dp.x)
  mort_rr_art_target <- array(NA, dim = c(3, 15, length(year.idx)), dimnames = list(transmission = c('0to6mo', '7to12mo', '12+mo'), age = 0:14, year = proj.years))
  mort_rr_art_target[1:2, 1:5,] <- rep(unlist(mort_rr_art[1,]), each = 10)
  mort_rr_art_target[3, 1:5,] <- rep(unlist(mort_rr_art[2,]), each = 5)
  mort_rr_art_target[1:2, 6:15,] <- rep(unlist(mort_rr_art[3,]), each = 20)
  mort_rr_art_target[3, 6:15,] <- rep(unlist(mort_rr_art[4,]), each = 10)
  v$mort_art_rr <- mort_rr_art_target

  art_dist_paed <- dp_read_art_dist(dp.x)
  v$init_art_dist <- art_dist_paed

  ##BF duration
  bf_duration <- input_breastfeeding_dur(dp.x)
  v$bf_duration_art <- bf_duration[,,2]
  v$bf_duration_no_art <- bf_duration[,,1]
  ##only keeping this for leapfrog
  v$bf_duration = bf_duration

  art_elig = dp_read_paed_art_eligibility(dp.x)
  v$paed_art_elig_age <- as.integer(art_elig$age_elig / 12) ##converts from months to years

  cd4_elig <- art_elig$cd4_elig[c(5:7,4),]
  ##Changing the input from CD4 count or percentages to ordinal categories
  ###Easier to do it here than in the leapfrog code
  paed_cd4_percent_intervals <- c(31, 30, 25, 20, 15, 10, 5)
  paed_cd4_number_intervals <- c(1001, 1000, 750, 500, 350, 200)
  paed_art_elig_cd4 <- array(data = NA, dim = c(15, length(year.idx)), dimnames = list(age = c(0:14), year = c(proj.years)))
  paed_art_elig_cd4[1,] <- findInterval(-unname(cd4_elig[1,]), -paed_cd4_percent_intervals)
  paed_art_elig_cd4[2:3,] <- findInterval(-unname(cd4_elig[2,]), -paed_cd4_percent_intervals)
  paed_art_elig_cd4[4:5,] <- findInterval(-unname(cd4_elig[3,]), -paed_cd4_percent_intervals)
  paed_art_elig_cd4[6:15,] <- rep(findInterval(-unname(cd4_elig[4,]), -paed_cd4_number_intervals), each = length(6:15))
  v$paed_art_elig_cd4 <- paed_art_elig_cd4


  v$paed_art_ltfu <- input_childart_ltfu(dp.x) / 100

  v$abortion <- input_abortion(dp.x)

  v$patients_reallocated <- input_mothers_reallocated(dp.x)


  ##extract needed outputs to just run paed model
  wlhiv_births <- dpsub(dp = dp.x, "<ChildNeedPMTCT MV>", 2, timedat.idx) %>% unname()
  names(wlhiv_births) <- proj.years
  rownames(wlhiv_births) <- NULL
  v$mat_hiv_births <- as.array(as.numeric(unlist(wlhiv_births)))
  v$mat_prev_input = rep(TRUE, length(year.idx))

  hivpop <- SpectrumUtils::dp.output.hivpop(dp.raw = dp.x, direction = 'long') %>% dplyr::rename(hivpop = Value)
  totpop <- SpectrumUtils::dp.output.bigpop(dp.raw = dp.x, direction = 'long') %>% dplyr::rename(totpop = Value)
  inc <- SpectrumUtils::dp.output.incident.hiv(dp.raw = dp.x, direction = 'long') %>% dplyr::rename(inc = Value)

  dt <- merge(hivpop, totpop, by = c('Sex', 'Age', 'Year'))
  dt <- dt %>% dplyr::mutate(hivnpop = totpop - hivpop) %>%
    dplyr::filter(Sex == 'Female' & Age %in% 15:49) %>%
    dplyr::select(Age, Year, hivnpop) %>%
    tidyr::pivot_wider(names_from = Year, values_from = hivnpop)

  hivnpop <- array(NA, dim = c(length(15:49), length(year.idx)), dimnames = list(age = 15:49, year = proj.years))
  for(i in 1:length(year.idx)){
    hivnpop[,i] <- dt[[(i+1)]]
  }

  inc <- inc %>%
    dplyr::filter(Sex == 'Female' & Age %in% 15:49) %>%
    dplyr::select(Age, Year, inc) %>%
    tidyr::pivot_wider(names_from = Year, values_from = inc)

  inc.array <- array(NA, dim = c(length(15:49), length(year.idx)), dimnames = list(age = 15:49, year = proj.years))
  for(i in 1:length(year.idx)){
    inc.array[,i] <- inc[[(i+1)]]
  }

  v$hivnpop_full <- hivnpop
  v$adult_female_infections_full <- inc.array

  h.fert.idx <- which((15L-1 + cumsum(params$hAG_SPAN_coarse)) %in% 15:49)
  fert_rat.h.ag <- findInterval(as.integer(rownames(hivnpop)),
                                15L + cumsum(v$hAG_SPAN_coarse[h.fert.idx]) -
                                  v$hAG_SPAN_coarse[h.fert.idx])
  coarse_age_groups <- cut(15:49,
                           breaks = 15L + cumsum(v$hAG_SPAN_coarse[h.fert.idx]) -
                             v$hAG_SPAN_coarse[h.fert.idx],
                           right = FALSE)
  v$hivnpop_coarse <- as.array(rowsum(hivnpop, group = coarse_age_groups))
  v$adult_female_infections_coarse <- as.array(rowsum(inc.array, group = coarse_age_groups))

  ##Make a coarse ASFR, needed to run WLHIV births at the coarse level
  fert_ages.idx <- 1L + cumsum(params$hAG_SPAN_coarse[h.fert.idx]) - params$hAG_SPAN_coarse[h.fert.idx]
  v$asfr_coarse <- as.array(rowsum(params$asfr, group = coarse_age_groups))

  total_births <- SpectrumUtils::dp.output.births(dp.raw = dp.x, direction = 'long')$Value %>% as.array()
  v$total_births <- total_births


  specres <- eppasm::read_hivproj_output(pjnz)
  newinf <- specres$newinf.f[4:10,] %>% colSums()
  newinf_rate <- newinf / colSums(specres$totpop.f[1:10,])

  fp1 <- eppasm::prepare_directincid(pjnz)
  fp1$tARTstart <- 61L
  #Need to run the adult model to pull out the propotion by CD4 category
  mod1 <- eppasm::simmod(fp1)

  #cd4
  wlhiv_cd4 <- array(as.numeric(unlist(dpsub(dp = dp.x, "<CD4Distribution15_49 MV2>", 19:25, timedat.idx))), dim = c(7,length(timedat.idx)))
  v$prop_gte350 <- colSums(wlhiv_cd4[1:2,]) / colSums(wlhiv_cd4)
  v$prop_lt200 <- colSums(wlhiv_cd4[5:7,]) / colSums(wlhiv_cd4)


  return(v)
}

df_as_numeric <- function(df) {
  vapply(df, function(x) as.numeric(x), numeric(nrow(df)))
}

df_as_logical <- function(df) {
  vapply(df, function(x) as.logical(as.numeric(x)), logical(nrow(df)))
}
