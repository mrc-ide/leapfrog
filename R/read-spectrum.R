
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
#' pjnz <- system.file("pjnz/bwa2021_v6.13.pjnz", package = "leapfrog")
#' dp <- read_dp(pjnz)
#' class(dp)
#'
#' @export
read_dp <- function(pjnz) {

  stopifnot(grepl("\\.(pjnz|zip)$", pjnz, ignore.case = TRUE))

  dpfile <- grep("\\.DP$", unzip(pjnz, list = TRUE)$Name, value = TRUE)
  stopifnot(length(dpfile) == 1)

  dp <- read.csv(unz(pjnz, dpfile), as.is = TRUE)
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

dpsub <- function(dp, tag, rows, cols, tagcol = 1) {

  stopifnot(inherits(dp, "spectrum_dp"))
  stopifnot(is.character(tag))
  all.equal(rows, as.integer(rows))
  all.equal(cols, as.integer(cols))

  dp[which(dp[, tagcol] == tag) + rows, cols]
}

get_dp_years <- function(dp) {
  yr_start <- as.integer(dpsub(dp, "<FirstYear MV2>", 2, 4))
  yr_end <- as.integer(dpsub(dp, "<FinalYear MV2>", 2, 4))

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
#' @examples
#'
#' pjnz <- system.file("pjnz/bwa2021_v6.13.pjnz", package = "leapfrog")
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
#'
#' @export
dp_read_anc_testing <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  anc_indicators <- c("anc_clients", "anc_tested", "anc_tested_pos", "anc_known_pos",
                      "anc_prevalence", "anc_retested", "anc_retested_pos",
                      "facility_births", "anc_known_neg")

  if (exists_dptag(dp, "<ANCTestingValues MV>")) {
    anc_testing <- dpsub(dp, "<ANCTestingValues MV>", 2:5, dpy$time_data_idx)
    anc_testing <- sapply(anc_testing, as.numeric)
    dimnames(anc_testing) <- list(indicator = anc_indicators[1:4], year = dpy$proj_years)
  } else if (exists_dptag(dp, "<ANCTestingValues MV2>")) {
    anc_testing <- dpsub(dp, "<ANCTestingValues MV2>", 2:5, dpy$time_data_idx)
    anc_testing <- sapply(anc_testing, as.numeric)
    dimnames(anc_testing) <- list(indicator = anc_indicators[1:4], year = dpy$proj_years)
  } else if (exists_dptag(dp, "<ANCTestingValues MV4>")) {
    anc_testing <- dpsub(dp, "<ANCTestingValues MV4>", 2:10, dpy$time_data_idx)
    anc_testing <- sapply(anc_testing, as.numeric)
    dimnames(anc_testing) <- list(indicator = anc_indicators, year = dpy$proj_years)
  } else {
    stop("ANC testing tag not recognized. Function probably needs update for this .DP file.")
  }

  anc_testing[anc_testing == -9999] <- NA_real_

  anc_testing
}

#' @rdname dp_read_anc_testing
#' @export
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
    pmtct_arv <- sapply(pmtct_arv, as.numeric)
    dimnames(pmtct_arv) <- list(indicator = pmtct_indicators, year = dpy$proj_years)
  } else {
    stop("PMTCT ARV tag not recognized. Function probably needs update for this .DP file.")
  }

  pmtct_arv[pmtct_arv == 0.0] <- NA_real_

  pmtct_arv
}


#' @rdname dp_read_anc_testing
#' @export
dp_read_pmtct_retained <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  indicator_names <- c("pmtct_retained_alreadyart", "pmtct_retained_newart")

  if (exists_dptag(dp, "<PercentARTDelivery MV>")) {
    pmtct_retained <- dpsub(dp, "<PercentARTDelivery MV>", 2:3, dpy$time_data_idx)
  } else {
    stop("PMTCT retained at delivery tag not recognized. Function probably needs update for this .DP file.")
  }

  pmtct_retained <- sapply(pmtct_retained, as.numeric)
  dimnames(pmtct_retained) <- list(indicator = indicator_names, year = dpy$proj_years)

  pmtct_retained
}

#' @rdname dp_read_anc_testing
#' @export
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

  pregtermabortion <- setNames(as.numeric(pregtermabortion), dpy$proj_years)
  pregtermabortion_ispercent <- setNames(as.logical(pregtermabortion_ispercent), dpy$proj_years)

  list(pregtermabortion = pregtermabortion,
       pregtermabortion_ispercent = pregtermabortion_ispercent)
}


#' @rdname dp_read_anc_testing
#' @export
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

  notbreastfeeding_percent <- sapply(notbreastfeeding_percent, as.numeric)
  notbreastfeeding_percent_noarv <- notbreastfeeding_percent[1:18, ]
  notbreastfeeding_percent_arv <- notbreastfeeding_percent[18 + 1:18, ]

  dn <- list(child_age_months = agemonth_cat, year = dpy$proj_years)
  dimnames(notbreastfeeding_percent_noarv) <- dn
  dimnames(notbreastfeeding_percent_arv) <- dn

  list(notbreastfeeding_percent_noarv = notbreastfeeding_percent_noarv,
       notbreastfeeding_percent_arv = notbreastfeeding_percent_arv)
}


#' @rdname dp_read_anc_testing
#' @export
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

  childart <- sapply(childart, as.numeric)
  dimnames(childart) <- list(indicator = indicator_names, year = dpy$proj_years)
  childart[childart == -9999] <- NA_real_

  if (exists_dptag(dp, "<ChildARTByAgeGroupPerNum MV2>")) {
    childart_ispercent <- dpsub(dp, "<ChildARTByAgeGroupPerNum MV2>", 2:6, dpy$time_data_idx)
  } else {
    stop("Child ART input tag not recognized. Function probably needs update for this .DP file.")
  }

  childart_ispercent <- sapply(childart_ispercent, as.logical)
  dimnames(childart_ispercent) <- list(indicator = indicator_names, year = dpy$proj_years)

  list(childart = childart,
       childart_ispercent = childart_ispercent)
}

#' @rdname dp_read_anc_testing
#' @export
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
#' @export
dp_read_art_dist <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  if (exists_dptag(dp, "<ChildARTDist MV>")) {
    child_mort_mult <- dpsub(dp, "<ChildARTDist MV>", 2:16, dpy$time_data_idx)
    child_mort_mult <- sapply(child_mort_mult, as.numeric)
    dimnames(child_mort_mult) <- list(age = 0:14, year = dpy$proj_years)
  } else {
    stop("Child ART distribution tag not recognized. Function probably needs update for this .DP file.")
  }

  child_mort_mult
}

#' @rdname dp_read_anc_testing
#' @export
dp_read_child_mort_mult <- function(dp) {

  dp <- get_dp_data(dp)
  dpy <- get_dp_years(dp)

  indicator_names <- c('Age <5, <12 months on ART',
                       'Age <5, 12+ months on ART',
                       'Age >=5, <12 months on ART',
                       'Age >=5, 12+ months on ART')

  if (exists_dptag(dp, "<ChildMortalityRates MV2>")) {
    child_mort_mult <- dpsub(dp, "<ChildMortalityRates MV2>", 2:5, dpy$time_data_idx)
    child_mort_mult <- sapply(child_mort_mult, as.numeric)
    dimnames(child_mort_mult) <- list(mrr = indicator_names, year = dpy$proj_years)
  } else {
    stop("Child ART multiplier tag not recognized. Function probably needs update for this .DP file.")
  }

  child_mort_mult
}
