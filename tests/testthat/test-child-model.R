test_that("Child model can be run for all years", {
  input <- setup_childmodel(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$parameters

  expect_silent(out <- run_model(demp, parameters, 1970:2030, NULL, run_child_model = TRUE))

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths", "p_hiv_pop",
      "p_hiv_pop_natural_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths", "hc1_hiv_pop", "hc2_hiv_pop",
      "hc1_art_pop", "hc2_art_pop",
      "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
      "hc1_art_aids_deaths", "hc2_art_aids_deaths", "hiv_births",
      "hc_art_init", "hc_art_need_init", "ctx_need",
      "ctx_mean"
    )
  )

  ## Nothing should ever be negative
  expect_true(all(out$hc1_hiv_pop[, , , , ] >= 0))
  expect_true(all(out$hc2_hiv_pop[, , , , ] >= 0))
  expect_true(all(out$hc1_art_pop[, , , , ] >= 0))
  expect_true(all(out$hc2_art_pop[, , , , ] >= 0))
  expect_true(all(out$hc1_noart_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc2_noart_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc1_art_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hc2_art_aids_deaths[, , , , ] >= 0))
  expect_true(all(out$hiv_births >= 0))
})

test_that("Infections among children align", {
  input <- setup_childmodel(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$parameters
  dp <- input$dp
  pjnz <- input$pjnz

  out <- run_model(demp, parameters, 1970:2030, NULL, run_child_model = TRUE)
  inf_spec <- SpectrumUtils::dp.output.incident.hiv(dp.raw = dp)
  inf_spec <- inf_spec %>%
    dplyr::filter(Sex == "Male+Female") %>%
    dplyr::filter(Age != "80+") %>%
    dplyr::mutate(Age = as.numeric(Age)) %>%
    dplyr::filter(Age < 15) %>%
    reshape2::melt(id.vars = c("Sex", "Age")) %>%
    dplyr::rename(Year = variable, Spec = value) %>%
    dplyr::mutate(Year = as.numeric(as.character(Year))) %>%
    dplyr::select(Age, Year, Spec) %>%
    dplyr::as_tibble()

  lfrog <- out$p_infections[1:15, , ] %>%
    reshape2::melt() %>%
    dplyr::rename(Age = Var1, Sex = Var2, Year = Var3, lfrog = value) %>%
    dplyr::mutate(Age = Age - 1, Year = Year + 1969) %>%
    dplyr::group_by(Age, Year) %>%
    dplyr::summarise(lfrog = sum(lfrog))

  dt <- dplyr::right_join(inf_spec, lfrog, by = c("Age", "Year"))
  dt <- dt %>%
    dplyr::mutate(diff = Spec - lfrog)

  expect_true(all(abs(dt$diff) < 5e-1))
})

test_that("CLHIV align", {
  input <- setup_childmodel(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$parameters
  dp <- input$dp
  pjnz <- input$pjnz

  out <- run_model(demp, parameters, 1970:2030, NULL, run_child_model = TRUE)

  spec_prev <- input$offtrt

  hc1 <- dplyr::right_join((reshape2::melt(out$hc1_hiv_pop)), data.frame(Var1 = 1:7, cd4_cat = c("gte30", "26-30", "21-25", "16-20", "11-15", "5-10", "lte5")), by = "Var1")
  hc2 <- dplyr::right_join(data.frame(reshape2::melt(out$hc2_hiv_pop)), data.frame(Var1 = 1:6, cd4_cat = c("gte1000", "750-999", "500-749", "350-499", "200-349", "lte200")), by = "Var1")
  hc1 <- dplyr::right_join(hc1, data.frame(Var2 = 1:4, transmission = c("perinatal", "bf0-6", "bf7-12", "bf12+")), by = "Var2")
  hc2 <- dplyr::right_join(hc2, data.frame(Var2 = 1:4, transmission = c("perinatal", "bf0-6", "bf7-12", "bf12+")), by = "Var2")
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
    dplyr::summarise(fr = sum(fr))
  dt <- dplyr::right_join(hc, spec_prev, by = c("sex", "age", "cd4_cat", "year", "transmission"))
  dt <- dt %>%
    dplyr::mutate(diff = pop - fr) %>%
    dplyr::filter(year < 2030)

  expect_true(all(abs(dt$diff) < 5e-1))
})

test_that("CLHIV on ART align", {
  input <- setup_childmodel(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$parameters
  dp <- input$dp
  pjnz <- input$pjnz

  out <- run_model(demp, parameters, 1970:2030, NULL, run_child_model = TRUE)

  spec_prev <- input$ontrt
  spec_prev <- spec_prev %>%
    dplyr::rename(time_art = transmission)
  hc1 <- dplyr::right_join(data.frame(reshape2::melt(out$hc1_art_pop)), data.frame(Var2 = 1:7, cd4_cat = c("gte30", "26-30", "21-25", "16-20", "11-15", "5-10", "lte5")), by = ("Var2"))
  hc1 <- dplyr::right_join(hc1, data.frame(Var1 = 1:3, time_art = c("ARTlte5mo", "ART6to12mo", "ARTgte12mo")), by = ("Var1"))
  hc2 <- dplyr::right_join(data.frame(reshape2::melt(out$hc2_art_pop)), data.frame(Var2 = 1:6, cd4_cat = c("gte1000", "750-999", "500-749", "350-499", "200-349", "lte200")), by = "Var2")
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
    dplyr::filter(year < 2030)

  expect_true(all(abs(dt$diff) < 5e-1))
})

test_that("HIV related deaths among CLHIV not on ART align", {
  input <- setup_childmodel(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$parameters
  dp <- input$dp
  pjnz <- input$pjnz
  aids_deathsnoart <- input$deaths_noart

  out <- run_model(demp, parameters, 1970:2030, NULL, run_child_model = TRUE)


  hc1 <- apply(out$hc1_noart_aids_deaths, c(3:5), sum)
  hc2 <- apply(out$hc2_noart_aids_deaths, c(3:5), sum)
  hc <- array(0, dim = c(15, 2, 61))
  hc[1:5, , ] <- hc1
  hc[6:15, , ] <- hc2
  dt <- dplyr::right_join(reshape2::melt(hc), reshape2::melt(aids_deathsnoart), by = c("Var1", "Var2", "Var3"))
  dt <- dt %>%
    dplyr::mutate(
      age = Var1 - 1,
      sex = dplyr::if_else(Var2 == 1, "Male", "Female"),
      year = Var3 + 1969,
      fr = value.x,
      spec = value.y
    ) %>%
    dplyr::select(age, sex, year, fr, spec)
  dt <- dt %>%
    dplyr::mutate(diff = spec - fr)

  expect_true(all(abs(dt$diff) < 5e-1))
})

test_that("HIV related deaths among CLHIV on ART align", {
  input <- setup_childmodel(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$parameters
  dp <- input$dp
  pjnz <- input$pjnz
  aids_deathsart <- input$deaths_art

  out <- run_model(demp, parameters, 1970:2030, NULL, run_child_model = TRUE)


  ## right now this is only working for the first year of ART, assuming its something with the timing on art
  hc1 <- apply(out$hc1_art_aids_deaths, c(3:5), sum)
  hc2 <- apply(out$hc2_art_aids_deaths, c(3:5), sum)
  hc <- array(0, dim = c(15, 2, 61))
  hc[1:5, , ] <- hc1
  hc[6:15, , ] <- hc2
  dt <- dplyr::right_join(reshape2::melt(hc), reshape2::melt(aids_deathsart), by = c("Var1", "Var2", "Var3"))
  dt <- dt %>%
    dplyr::mutate(
      age = Var1 - 1,
      sex = dplyr::if_else(Var2 == 1, "Male", "Female"),
      year = Var3 + 1969,
      fr = value.x,
      spec = value.y
    ) %>%
    dplyr::select(age, sex, year, fr, spec)
  dt <- dt %>%
    dplyr::mutate(diff = spec - fr) %>%
    dplyr::filter(year < 2030)

  expect_true(all(abs(dt$diff) < 5e-1))
})

