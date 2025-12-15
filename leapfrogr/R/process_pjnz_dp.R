process_pjnz_dp <- function(dat, pars, dim_vars) {
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
  dimnames(pars$asfr) <- list(age = 15:49, year = dim_vars$years)

  pars$migr_age_dist <- pars$migr_age_dist / 100
  migr_age_dist_sum <- colSums(pars$migr_age_dist)
  migr_age_dist_sum[migr_age_dist_sum == 0] <- 1
  pars$migr_age_dist <- sweep(pars$migr_age_dist, 2:3, migr_age_dist_sum, "/")

  netmigr_5year <- sweep(pars$migr_age_dist, 2:3, pars$migr_rate, "*")
  netmigr <- array(
    dim = c(length(dim_vars$a_80plus_no_80), length(dim_vars$g), length(dim_vars$years)),
    dimnames = list(
      age = dim_vars$a_80plus_no_80,
      sex = dim_vars$g,
      year = dim_vars$years
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

  pars
}
