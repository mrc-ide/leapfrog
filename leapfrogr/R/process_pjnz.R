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
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "leapfrog"
#' )
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

#' Prepare inputs from Spectrum PJNZ
#'
#' @param pjnz path to PJNZ file
#' @param use_coarse_age_groups use the coarse age stratification
#' @param bypass_adult produce parameters that will bypass the adult
#'  model when running the child model variant
#'
#' @return list of input parameters
#'
#' @examples
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "leapfrog")
#' parameters <- process_pjnz(pjnz)
#'
#' @export
process_pjnz <- function(pjnz, use_coarse_age_groups = FALSE, bypass_adult = FALSE) {
  dp <- read_dp(pjnz)
  dat <- parse_dp(dp)
  dim_vars <- dat$dim_vars

  pars <- lapply(dat$data, function(x) if (is.null(x)) NULL else x$data)
  names(pars) <- names(dat$data)

  pars <- process_pjnz_dp(dat, pars, dim_vars)
  pars <- process_pjnz_ha(dat, pars, dim_vars, dp, use_coarse_age_groups)
  pars <- process_pjnz_hc(dat, pars, dim_vars, use_coarse_age_groups, bypass_adult)
  pars
}
