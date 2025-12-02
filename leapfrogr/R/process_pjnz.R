process_pjnz <- function(pjnz, use_coarse_age_groups = FALSE, bypass_adult = FALSE) {
  dp <- read_dp(pjnz)
  dat <- parse_dp(dp)
  dim_vars <- dat$dim_vars

  pars <- lapply(dat$data, function(x) if (is.null(x)) NULL else x$data)
  names(pars) <- names(dat$data)

  pars <- process_pjnz_dp(dat, pars, dim_vars)
  pars <- process_pjnz_ha(dat, pars, dim_vars, use_coarse_age_groups, pjnz = pjnz)
  pars <- process_pjnz_hc(dat, pars, dim_vars, use_coarse_age_groups, bypass_adult)
  pars
}
