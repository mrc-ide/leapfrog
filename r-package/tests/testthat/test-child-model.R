test_that("Child model can be run for all years", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  expect_silent(out <- run_model(parameters, "ChildModel", 1970:2030))

  expect_setequal(
    names(out),
    c(
      "p_totpop", "births", "p_deaths_background_totpop", "p_hivpop",
      "p_deaths_background_hivpop", "h_hivpop", "h_artpop",
      "h_hiv_deaths_no_art", "h_deaths_excess_nonaids_no_art",
      "h_deaths_excess_nonaids_on_art", "p_deaths_excess_nonaids",
      "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths", "p_net_migration_hivpop",
      "hiv_births_by_mat_age", "hiv_births",
      "hc1_hivpop", "hc2_hivpop",
      "hc1_artpop", "hc2_artpop",
      "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
      "hc1_art_aids_deaths", "hc2_art_aids_deaths",
      "hc_art_init", "hc_art_need_init", "ctx_need", "infection_by_type",
      "mtct_by_source_tr", "mtct_by_source_women",
      "mtct_by_source_hc_infections", "pmtct_coverage_at_delivery")
  )

  ## Nothing should ever be negative
  expect_true(all(out$hiv_births >= 0))
  expect_true(all(out$hc1_hivpop[, , , , ] >= 0))
  expect_true(all(out$hc2_hivpop[, , , , ] >= 0))
  expect_true(all(out$hc1_artpop[, , , , ] >= 0))
  expect_true(all(out$hc2_artpop[, , , , ] >= 0))
  expect_true(all(out$hc1_noart_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc2_noart_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc1_art_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc2_art_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc1_hivpop_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc2_hivpop_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc1_artpop_strat[, , , , , , , ] >= 0))
  expect_true(all(out$hc2_artpop_strat[, , , , , , , ] >= 0))
  expect_true(all(out$hc1_noart_aids_deaths_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc2_noart_aids_deaths_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc1_art_aids_deaths_strat[, , , , , , , ] >= 0))
  expect_true(all(out$hc2_art_aids_deaths_strat[, , , , , , , ] >= 0))
  expect_true(all(out$pmtct_coverage_at_delivery >= 0))
})

test_that("Coarse child model can be run for all years", {
  parameters <- read_parameters(test_path("testdata/child_parms_coarse.h5"))
  parameters$mat_prev_input[] <- as.integer(0)
  out_coarse <- run_model(parameters, "CoarseChildModel", 1970:2030)

  expect_setequal(
    names(out_coarse),
    c(
      "p_totpop", "births", "p_deaths_background_totpop", "p_hivpop",
      "p_deaths_background_hivpop", "h_hivpop", "h_artpop",
      "h_hiv_deaths_no_art", "h_deaths_excess_nonaids_no_art",
      "h_deaths_excess_nonaids_on_art", "p_deaths_excess_nonaids",
      "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths", "p_net_migration_hivpop",
      "hiv_births_by_mat_age", "hiv_births",
      "hc1_hivpop", "hc2_hivpop",
      "hc1_artpop", "hc2_artpop",
      "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
      "hc1_art_aids_deaths", "hc2_art_aids_deaths",
      "hc_art_init", "hc_art_need_init", "ctx_need", "infection_by_type",
      "mtct_by_source_tr", "mtct_by_source_women",
      "mtct_by_source_hc_infections", "pmtct_coverage_at_delivery")
  )

  ## Nothing should ever be negative
  expect_true(all(out_coarse$hiv_births >= 0))
  expect_true(all(out_coarse$hiv_births_by_mat_age >= 0))
  expect_true(all(out_coarse$hc1_hivpop[, , , , ] >= 0))
  expect_true(all(out_coarse$hc2_hivpop[, , , , ] >= 0))
  expect_true(all(out_coarse$hc1_artpop[, , , , ] >= 0))
  expect_true(all(out_coarse$hc2_artpop[, , , , ] >= 0))
  expect_true(all(out_coarse$hc1_noart_aids_deaths[, , , , ] >= 0))
  expect_true(all(out_coarse$hc2_noart_aids_deaths[, , , , ] >= 0))
  expect_true(all(out_coarse$hc1_art_aids_deaths[, , , , ] >= 0))
  expect_true(all(out_coarse$hc2_art_aids_deaths[, , , , ] >= 0))
  expect_true(all(out_coarse$hc1_hivpop_strat[, , , , , , ] >= 0))
  expect_true(all(out_coarse$hc2_hivpop_strat[, , , , , , ] >= 0))
  expect_true(all(out_coarse$hc1_artpop_strat[, , , , , , , ] >= 0))
  expect_true(all(out_coarse$hc2_artpop_strat[, , , , , , , ] >= 0))
  expect_true(all(out_coarse$hc1_noart_aids_deaths_strat[, , , , , , ] >= 0))
  expect_true(all(out_coarse$hc2_noart_aids_deaths_strat[, , , , , , ] >= 0))
  expect_true(all(out_coarse$hc1_art_aids_deaths_strat[, , , , , , , ] >= 0))
  expect_true(all(out_coarse$hc2_art_aids_deaths_strat[, , , , , , , ] >= 0))
  expect_true(all(out_coarse$pmtct_coverage_at_delivery >= 0))
})

test_that("Model outputs are consistent", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  out <- run_model(parameters, "ChildModel", 1970:2030)

  ###############################
  ##HIV births by maternal age matches total hiv births
  ###############################
  if(!all(parameters$mat_prev_input)){
    expect_true(all(abs(colSums(out$hiv_births_by_mat_age) - out$hiv_births) < 1e-5))
  }

  ###############################
  ##Infections stratified by infection type and population infections should be the same
  ###############################
  strat <- apply(out$infection_by_type, c(2,3,4), sum)
  pop <- out$p_infections[1:5,,]
  expect_equal(strat, pop)

  ###############################
  ##Infections stratified by stacked bar categories should be the same
  ###############################
  strat <- apply(out$infection_by_type, c(1,4), sum)
  stacked_bar <- apply(out$mtct_by_source_hc_infections, c(2,3), sum)
  #perinatal infections
  expect_equal(strat[1,], stacked_bar[1,])
  #bf infections
  expect_equal(apply(strat[2:4,], 2, 'sum'), apply(out$mtct_by_source_hc_infections[,2:5,], 3, 'sum'))

  ###############################
  ##Stratified deaths and population deaths should be the same
  ###############################
  ##p_hiv_deaths & hc1_art_aids_deaths, hc1_noart_aids_deaths, hc2_art_aids_deaths, hc2_noart_aids_deaths
  hc1_hiv <- apply(out$hc1_noart_aids_deaths, c(3,4,5), sum)
  hc1_art <- apply(out$hc1_art_aids_deaths, c(3,4,5), sum)
  hc1 <- hc1_hiv + hc1_art
  dimnames(hc1) <- list(age = 0:4, sex = c('male','female'), year = 1970:2030)
  p_hiv_deaths <- out$p_hiv_deaths[1:5, , ]
  dimnames(p_hiv_deaths) <- list(age = 0:4, sex = c('male','female'), year = 1970:2030)

  hc1_df <- as.data.frame(as.table(hc1))
  out_df <- as.data.frame(as.table(p_hiv_deaths[1:5, , ]))
  colnames(hc1_df) <- c("Var1", "Var2", "Var3", "strat")
  colnames(out_df) <- c("Var1", "Var2", "Var3", "pop")
  c1 <- hc1_df %>%
    dplyr::inner_join(out_df, by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(diff = strat - pop)
  expect_true(all(abs(c1$diff) < 1e-5))

  hc2_hiv <- apply(out$hc2_noart_aids_deaths, c(3,4,5), sum)
  hc2_art <- apply(out$hc2_art_aids_deaths, c(3,4,5), sum)
  hc2 <- hc2_hiv + hc2_art
  dimnames(hc2) <- list(age = 5:14, sex = c('male','female'), year = 1970:2030)
  p_hiv_deaths <- out$p_hiv_deaths[6:15, , ]
  dimnames(p_hiv_deaths) <- list(age = 5:14, sex = c('male','female'), year = 1970:2030)
  hc2_df <- as.data.frame(as.table(hc2))
  out_df <- as.data.frame(as.table(p_hiv_deaths))
  colnames(hc2_df) <- c("Var1", "Var2", "Var3", "strat")
  colnames(out_df) <- c("Var1", "Var2", "Var3", "pop")
  c2 <- hc2_df %>%
    dplyr::inner_join(out_df, by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(diff = strat - pop)
  expect_true(all(abs(c2$diff) < 1e-5))


  ###############################
  ##Stratified hiv pop and population hiv pop should be the same
  ###############################
  ##p_hivpop & hc1_hivpop, hc1_artpop, hc2_hivpop, hc2_hivpop
  hc1_hiv <- apply(out$hc1_hivpop, c(3,4,5), sum)
  hc1_art <- apply(out$hc1_artpop, c(3,4,5), sum)
  hc1 <- hc1_hiv + hc1_art
  dimnames(hc1) <- list(age = 0:4, sex = c('male','female'), year = 1970:2030)
  p_hiv <- out$p_hivpop[1:5, , ]
  dimnames(p_hiv) <- list(age = 0:4, sex = c('male','female'), year = 1970:2030)
  hc1_df <- as.data.frame(as.table(hc1))
  out_df <- as.data.frame(as.table(p_hiv))
  colnames(hc1_df) <- c("Var1", "Var2", "Var3", "strat")
  colnames(out_df) <- c("Var1", "Var2", "Var3", "pop")
  c1 <- hc1_df %>%
    dplyr::inner_join(out_df, by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(diff = strat - pop)
  expect_true(all(abs(c1$diff) < 1e-5))

  hc2_hiv <- apply(out$hc2_hivpop, c(3,4,5), sum)
  hc2_art <- apply(out$hc2_artpop, c(3,4,5), sum)
  hc2 <- hc2_hiv + hc2_art
  dimnames(hc2) <- list(age = 5:14, sex = c('male','female'), year = 1970:2030)
  p_hiv <- out$p_hivpop[6:15, , ]
  dimnames(p_hiv) <- list(age = 5:14, sex = c('male','female'), year = 1970:2030)
  hc2_df <- as.data.frame(as.table(hc2))
  out_df <- as.data.frame(as.table(p_hiv))
  colnames(hc2_df) <- c("Var1", "Var2", "Var3", "strat")
  colnames(out_df) <- c("Var1", "Var2", "Var3", "pop")
  c2 <- hc2_df %>%
    dplyr::inner_join(out_df, by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(diff = strat - pop)
  expect_true(all(abs(c2$diff) < 1e-5))
})

test_that("Female 15-49y pop aligns", {
  testthat::skip("Skipping this test because the adult populations currently do not align")
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)
  spec <- SpectrumUtils::dp.output.hivpop(dp, direction = 'long')
  spec <- spec %>%
    dplyr::filter(Age %in% 15:49 & Sex == 'Female') %>%
    dplyr::group_by(Year)

  lfrog <- out$p_hivpop[16:50,2,]
  dimnames(lfrog) <- list(Age = 15:49, Year = 1970:2030)
  lfrog <- as.data.frame(as.table(lfrog)) %>%
     dplyr::mutate(lfrog = Freq)
  lfrog$Year <- as.integer(as.character(lfrog$Year))


  dt <- dplyr::right_join(spec, lfrog, by = c("Year", "Age"))
  dt <- dt %>%
    dplyr::mutate(diff = Value - lfrog)

  expect_true(all(abs(dt$diff) < 1e-3))
})

test_that("Mothers that need ptmct align", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)
  spec <- SpectrumUtils::dp.output.pmtct.need(dp, direction = 'long')

  lfrog <- data.frame(lfrog = out$hiv_births, Year = 1970:2030)

  dt <- dplyr::right_join(spec, lfrog, by = c("Year"))
  dt <- dt %>%
    dplyr::mutate(diff = Value - lfrog)

  expect_true(all(abs(dt$diff) < 1e-3))
})

test_that("Children in need of cotrim aligns", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  out <- run_model(parameters, "ChildModel", 1970:2030)

  spec <- utils$ctx_need
  dt <- data.frame(year = 1970:2030,
                   spec = as.numeric(unlist(spec)),
                   lfrog = out$ctx_need)
  dt <- dt %>%
    dplyr::mutate(diff = spec - lfrog)

  expect_true(all(abs(dt$diff) < 5e-1))
})

test_that("Infections among children align", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)
  system.time(out <- run_model(parameters, "ChildModel", 1970:2030))
  inf_spec <- SpectrumUtils::dp.output.incident.hiv(dp.raw = dp)
  inf_spec <- inf_spec %>%
    dplyr::filter(Sex != "Male+Female") %>%
    dplyr::filter(Age != "80+") %>%
    dplyr::mutate(Age = as.numeric(Age)) %>%
    dplyr::filter(Age < 15) %>%
    reshape2::melt(id.vars = c("Sex", "Age")) %>%
    dplyr::rename(Year = variable, Spec = value) %>%
    dplyr::mutate(Year = as.numeric(as.character(Year))) %>%
    dplyr::select(Age, Year, Sex, Spec) %>%
    dplyr::as_tibble()

  lfrog <- out$p_infections[1:15, , ] %>%
    reshape2::melt() %>%
    dplyr::rename(Age = Var1, Sex = Var2, Year = Var3, lfrog = value) %>%
    dplyr::mutate(Age = Age - 1, Year = Year + 1969, Sex = ifelse(Sex == 1, 'Male', 'Female')) %>%
    dplyr::group_by(Age, Year, Sex) %>%
    dplyr::summarise(lfrog = sum(lfrog), .groups = c("keep"))

  dt <- dplyr::right_join(inf_spec, lfrog, by = c("Age", "Year", "Sex"))
  dt <- dt %>%
    dplyr::mutate(diff = Spec - lfrog) %>%
    dplyr::filter(Age < 4 & Year < 2030)

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("CLHIV align", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)

  spec_prev <- utils$offtrt

  hc1 <- dplyr::right_join((reshape2::melt(out$hc1_hivpop)), data.frame(Var1 = 1:7, cd4_cat = c("gte30", "26-30", "21-25", "16-20", "11-15", "5-10", "lte5")), by = "Var1")
  hc2 <- dplyr::right_join(data.frame(reshape2::melt(out$hc2_hivpop)), data.frame(Var1 = 1:6, cd4_cat = c("gte1000", "750-999", "500-749", "350-499", "200-349", "lte200")), by = "Var1")
  hc1 <- dplyr::right_join(hc1, data.frame(Var2 = 1:5, transmission = c("perinatal", "bf0-6", "bf7-12", "bf12+", "bf12+")), by = "Var2")
  hc2 <- dplyr::right_join(hc2, data.frame(Var2 = 1:5, transmission = c("perinatal", "bf0-6", "bf7-12", "bf12+", "bf12+")), by = "Var2")
  hc1 <- hc1 %>%
    dplyr::mutate(
      age = Var3 - 1,
      sex = dplyr::if_else(Var4 == 1, "Male", "Female"),
      year = Var5 + 1969,
      fr = value
    ) %>%
    dplyr::select(age, cd4_cat, transmission, sex, year, fr)
  hc2 <- hc2 %>%
    dplyr::mutate(
      age = Var3 + 4,
      sex = dplyr::if_else(Var4 == 1, "Male", "Female"),
      year = Var5 + 1969,
      fr = value
    ) %>%
    dplyr::select(age, cd4_cat, transmission, sex, year, fr)
  hc <- rbind(hc1, hc2)
  hc <- hc %>%
    dplyr::mutate(
      age = as.integer(age),
      cd4_cat = as.factor(cd4_cat),
      sex = as.factor(sex),
      year = as.integer(year),
      transmission = as.factor(transmission)
    ) %>%
    dplyr::group_by(age, cd4_cat, sex, year, transmission) %>%
    dplyr::summarise(fr = sum(fr), .groups = c("keep"))
  dt <- dplyr::right_join(hc, spec_prev, by = c("sex", "age", "cd4_cat", "year", "transmission"))
  dt <- dt %>%
    dplyr::mutate(diff = pop - fr) %>%
    dplyr::filter(year < 2030) %>%
    dplyr::filter(age < 15)

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("CLHIV on ART align", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)

  spec_prev <- utils$ontrt
  spec_prev <- spec_prev %>%
    dplyr::rename(time_art = transmission)
  hc1 <- dplyr::right_join(data.frame(reshape2::melt(out$hc1_artpop)), data.frame(Var2 = 1:7, cd4_cat = c("gte30", "26-30", "21-25", "16-20", "11-15", "5-10", "lte5")), by = ("Var2"))
  hc1 <- dplyr::right_join(hc1, data.frame(Var1 = 1:3, time_art = c("ARTlte5mo", "ART6to12mo", "ARTgte12mo")), by = ("Var1"))
  hc2 <- dplyr::right_join(data.frame(reshape2::melt(out$hc2_artpop)), data.frame(Var2 = 1:6, cd4_cat = c("gte1000", "750-999", "500-749", "350-499", "200-349", "lte200")), by = "Var2")
  hc2 <- dplyr::right_join(hc2, data.frame(Var1 = 1:3, time_art = c("ARTlte5mo", "ART6to12mo", "ARTgte12mo")), by = ("Var1"))
  hc1 <- hc1 %>%
    dplyr::mutate(
      age = Var3 - 1,
      sex = dplyr::if_else(Var4 == 1, "Male", "Female"),
      year = Var5 + 1969,
      fr = value
    ) %>%
    dplyr::select(age, cd4_cat, time_art, sex, year, fr)
  hc2 <- hc2 %>%
    dplyr::mutate(
      age = Var3 + 4,
      sex = dplyr::if_else(Var4 == 1, "Male", "Female"),
      year = Var5 + 1969,
      fr = value
    ) %>%
    dplyr::select(age, cd4_cat, time_art, sex, year, fr)
  hc <- rbind(hc1, hc2)
  hc <- hc %>%
    dplyr::mutate(
      age = as.integer(age),
      cd4_cat = as.factor(cd4_cat),
      time_art = as.factor(time_art),
      sex = as.factor(sex),
      year = as.integer(year)
    ) %>%
    dplyr::group_by(age, cd4_cat, time_art, sex, year)
  dt <- dplyr::right_join(hc, spec_prev, by = c("sex", "age", "cd4_cat", "year", "time_art"))
  dt <- dt %>%
    dplyr::mutate(diff = pop - fr) %>%
    dplyr::filter(year < 2030 & age < 15)

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("HIV related deaths among CLHIV not on ART align", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz
  aids_deathsnoart <- utils$deaths_noart

  out <- run_model(parameters, "ChildModel", 1970:2030)

  hc1 <- apply(out$hc1_noart_aids_deaths, c(3:5), sum)
  hc2 <- apply(out$hc2_noart_aids_deaths, c(3:5), sum)
  hc <- array(0, dim = c(15, 2, 61))
  hc[1:5, , ] <- hc1
  hc[6:15, , ] <- hc2
  dimnames(hc) <- dimnames(aids_deathsnoart)
  dt <- dplyr::right_join(reshape2::melt(hc), reshape2::melt(aids_deathsnoart), by = c("age", "sex", "years"))
  dt <- dt %>%
    dplyr::mutate(
      year = years,
      fr = value.x,
      spec = value.y
    ) %>%
    dplyr::select(age, sex, year, fr, spec)
  dt <- dt %>%
    dplyr::mutate(diff = spec - fr)

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("HIV related deaths among CLHIV on ART align", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz
  aids_deathsart <- utils$deaths_art

  out <- run_model(parameters, "ChildModel", 1970:2030)

  hc1 <- apply(out$hc1_art_aids_deaths, c(3:5), sum)
  hc2 <- apply(out$hc2_art_aids_deaths, c(3:5), sum)
  hc <- array(0, dim = c(15, 2, 61))
  hc[1:5, , ] <- hc1
  hc[6:15, , ] <- hc2
  dimnames(hc) <- dimnames(aids_deathsart)
  dt <- dplyr::right_join(reshape2::melt(hc), reshape2::melt(aids_deathsart), by = c("age", "sex", "years"))
  dt <- dt %>%
    dplyr::mutate(
      year = years,
      fr = value.x,
      spec = value.y
    ) %>%
    dplyr::select(age, sex, year, fr, spec)
  dt <- dt %>%
    dplyr::mutate(diff = spec - fr) %>%
    dplyr::filter(year < 2030)

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("Child model agrees when run through all years vs two parts vs single year runs", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))

  # All years
  out_all_years <- run_model(parameters, "ChildModel", 1970:2030)

  # In two parts
  out_first_half_years <- run_model(parameters, "ChildModel", 1970:2000)
  out_second_half_years <- run_model_from_state(parameters, "ChildModel", get_time_slice(out_first_half_years, 31), 2000, 2001:2030)

  expect_equal(out_all_years, concat_on_time_dim(out_first_half_years, out_second_half_years))

  # Single years
  out_single_year <- get_time_slice(run_model(parameters, "ChildModel", 1970), 1)
  for (year in 1971:2030) {
    out_single_year <- run_model_single_year(parameters, "ChildModel", out_single_year, year - 1)
    expect_equal(out_single_year, get_time_slice(out_all_years, year - 1970 + 1))
  }
})

test_that("error thrown if trying to output from years before simulation start year", {
  parameters <- read_parameters(test_path("testdata/child_parms_full.h5"))

  out_first_half_years <- run_model(parameters, "ChildModel", 1970:2000)
  expect_error(
    run_model_from_state(parameters, "ChildModel", get_time_slice(out_first_half_years, 31), 2000, 2000:2030),
    "Cannot output year '2000'. Output years must be later than simulation start year '2000'.",
    fixed = TRUE
  )
})

