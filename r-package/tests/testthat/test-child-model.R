test_that("Child model can be run for all years", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))

  expect_silent(out <- run_model(parameters, "ChildModel", 1970:2030))

  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_background_deaths", "p_hiv_pop",
      "p_hiv_pop_background_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths", "hc1_hiv_pop", "hc2_hiv_pop",
      "hc1_art_pop", "hc2_art_pop",
      "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
      "hc1_art_aids_deaths", "hc2_art_aids_deaths", "hiv_births",
      "hc_art_init", "hc_art_need_init", "ctx_need", "infection_by_type",
      "mtct_by_source_tr", "mtct_by_source_women", "mtct_by_source_hc_infections")
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

  expect_true(all(out$hc1_hiv_pop_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc2_hiv_pop_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc1_art_pop_strat[, , , , , , , ] >= 0))
  expect_true(all(out$hc2_art_pop_strat[, , , , , , , ] >= 0))
  expect_true(all(out$hc1_noart_aids_deaths_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc2_noart_aids_deaths_strat[, , , , , , ] >= 0))
  expect_true(all(out$hc1_art_aids_deaths_strat[, , , , , , , ] >= 0))
  expect_true(all(out$hc2_art_aids_deaths_strat[, , , , , , , ] >= 0))

})

test_that("Model outputs are consistent", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
  out <- run_model(parameters, "ChildModel", 1970:2030)

  ###############################
  ##Infections stratified by infection type and population infections should be the same
  ###############################
  strat <- apply(out$infection_by_type, c(2,3,4), sum)
  pop <- out$p_infections[1:5,,]
  expect_equal(strat, pop)

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
  ##p_hiv_pop & hc1_hiv_pop, hc1_art_pop, hc2_hiv_pop, hc2_hiv_pop
  hc1_hiv <- apply(out$hc1_hiv_pop, c(3,4,5), sum)
  hc1_art <- apply(out$hc1_art_pop, c(3,4,5), sum)
  hc1 <- hc1_hiv + hc1_art
  dimnames(hc1) <- list(age = 0:4, sex = c('male','female'), year = 1970:2030)
  p_hiv <- out$p_hiv_pop[1:5, , ]
  dimnames(p_hiv) <- list(age = 0:4, sex = c('male','female'), year = 1970:2030)
  hc1_df <- as.data.frame(as.table(hc1))
  out_df <- as.data.frame(as.table(p_hiv))
  colnames(hc1_df) <- c("Var1", "Var2", "Var3", "strat")
  colnames(out_df) <- c("Var1", "Var2", "Var3", "pop")
  c1 <- hc1_df %>%
    dplyr::inner_join(out_df, by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(diff = strat - pop)
  expect_true(all(abs(c1$diff) < 1e-5))

  hc2_hiv <- apply(out$hc2_hiv_pop, c(3,4,5), sum)
  hc2_art <- apply(out$hc2_art_pop, c(3,4,5), sum)
  hc2 <- hc2_hiv + hc2_art
  dimnames(hc2) <- list(age = 5:14, sex = c('male','female'), year = 1970:2030)
  p_hiv <- out$p_hiv_pop[6:15, , ]
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
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)
  spec <- SpectrumUtils::dp.output.hivpop(dp, direction = 'long')
  spec <- spec %>%
    dplyr::filter(Age %in% 15:49 & Sex == 'Female') %>%
    dplyr::group_by(Year)

  lfrog <- out$p_hiv_pop[16:50,2,]
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
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
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
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
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

test_that("Stacked bar outputs align",{
  warning('MTCT by source output from Spectrum must be manually extracted when the Spectrum file is changed')
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
  out <- run_model(parameters, "ChildModel", 1970:2030)
  spec_output <- read.csv(test_path("testdata/mtct_by_source_MANUALLY_GENERATED.csv"))

  ###Compare number of women
  dimnames(out$mtct_by_source_women) <- list(Source = c(rep('Started ART during the pregnancy', 4),
                                                        'Started ART before the pregnancy',
                                                        'Started ART during the pregnancy',
                                                        'Started ART late in the pregnancy',
                                                        'Did not receive ART',
                                                        'Mother infected during pregnancy',
                                                        'Started ART before pregnancy then dropped off',
                                                        'Started ART during pregnancy then dropped off'),
                                             Year = 1970:2030)
  frogger_output_women <- reshape2::melt(out$mtct_by_source_women) %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(Source, Year) %>%
    dplyr::summarise(frogger = sum(value), .groups = "drop")
  spec_output_women <- spec_output %>%
    dplyr::filter(Indicator == "Number of pregnant women") %>%
    dplyr::transmute(Source, Year, spectrum = Value)
  women <- merge(frogger_output_women, spec_output_women)
  women <- dplyr::mutate(women, diff = frogger - spectrum)
  expect_true(all(abs(women$diff) < 5e-4))

  ###Compare transmission rate
  dimnames(out$mtct_by_source_tr) <- list(Source = c(rep('Started ART during the pregnancy', 4),
                                                     'Started ART before the pregnancy',
                                                     'Started ART during the pregnancy',
                                                     'Started ART late in the pregnancy',
                                                     'Did not receive ART',
                                                     'Mother infected',
                                                     'Started ART before pregnancy then dropped off',
                                                     'Started ART during pregnancy then dropped off'),
                                          trans_type = c('Perinatal', 'Breastfeeding'),
                                          Year = 1970:2030)
  frogger_output_tr <- out$mtct_by_source_tr %>%
    as.table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::rename(value = Freq) %>%
    dplyr::mutate(
      Source = dplyr::if_else(
        Source %in% c("Started ART during pregnancy then dropped off",
                      "Started ART before pregnancy then dropped off"),
        "Dropped off ART",
        Source
      )
    ) %>%
    dplyr::group_by(Source, trans_type, Year) %>%
    dplyr::summarise(frogger = sum(value), .groups = "drop")
  spec_output_tr <- spec_output %>% dplyr::filter(Indicator == "Mother-to-child transmission rate")
  spec_output_tr <- spec_output_tr %>% dplyr::mutate(spectrum = Value) %>%
    dplyr::mutate(trans_type = ifelse(grepl("breastfeeding", Source), "Breastfeeding", "Perinatal")) %>%
    dplyr::select(Source, trans_type, Year, spectrum) %>%
    dplyr::mutate(Source = gsub(pattern = ': child infected during pregnancy', replacement = '', x = spec_output_tr$Source)) %>%
    dplyr::mutate(Source = gsub(pattern = ': child infected during breastfeeding', replacement = '', x = spec_output_tr$Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Started ART before', Source), 'Started ART before the pregnancy', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Started ART during', Source), 'Started ART during the pregnancy', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Started ART late', Source), 'Started ART late in the pregnancy', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Dropped off ART during pregnancy', Source), 'Dropped off ART', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Mother infected', Source), 'Mother infected', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Did not receive ART', Source), 'Did not receive ART', Source)) %>%
    dplyr::mutate(spectrum  = spectrum / 100)
  tr <- merge(frogger_output_tr, spec_output_tr, by = c('Source', 'trans_type', 'Year'))
  tr <- dplyr::mutate(tr, diff = frogger - spectrum)
  expect_true(all(abs(tr$diff) < 5e-3))

  ###Ensure that the total frogger transmission rate matches with the stacked bar
  tr <- tr %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(Year, trans_type) %>%
    dplyr::summarise(frogger = sum(frogger), .groups = "drop") %>%
    dplyr::mutate(Year = as.integer(Year))
  dimnames(out$infection_by_type) <- list(
    trans_type = c('Perinatal', rep('Breastfeeding', 3)),
    age = 0:4,
    sex = c('Male', 'Female'),
    Year = 1970:2030
  )
  inf <- reshape2::melt(out$infection_by_type) %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(trans_type, Year) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    dplyr::left_join(
      tibble::tibble(hiv_births = out$hiv_births, Year = 1970:2030),
      by = "Year"
    ) %>%
    dplyr::mutate(tr = value / hiv_births)
  tr <- tr %>%
    dplyr::left_join(inf, by = c("Year", "trans_type")) %>%
    dplyr::mutate(diff = frogger - tr)
  expect_true(all(abs(tr$diff) < 5e-5))

  dimnames(out$mtct_by_source_hc_infections) <- list(Source = c(rep('Started ART during the pregnancy', 4),
                                                                'Started ART before the pregnancy',
                                                                'Started ART during the pregnancy',
                                                                'Started ART late in the pregnancy',
                                                                'Did not receive ART',
                                                                'Mother infected',
                                                                'Started ART before pregnancy then dropped off',
                                                                'Started ART during pregnancy then dropped off'),
                                          trans_type = c('Perinatal', 'Breastfeeding'),
                                          Year = 1970:2030)
  frogger_output_infections <- out$mtct_by_source_hc_infections %>%
    as.table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::rename(value = Freq) %>%
    dplyr::mutate(
      Source = dplyr::if_else(
        Source %in% c("Started ART during pregnancy then dropped off",
                      "Started ART before pregnancy then dropped off"),
        "Dropped off ART",
        Source
      )
    ) %>%
    dplyr::group_by(Source, trans_type, Year) %>%
    dplyr::summarise(frogger = round(sum(value)), .groups = "drop")
  spec_output_infections <- spec_output %>% dplyr::filter(Indicator == "New child infections")
  spec_output_infections <- spec_output_infections %>% dplyr::mutate(spectrum = Value) %>%
    dplyr::mutate(trans_type = ifelse(grepl("breastfeeding", Source), "Breastfeeding", "Perinatal")) %>%
    dplyr::select(Source, trans_type, Year, spectrum) %>%
    dplyr::mutate(Source = gsub(pattern = ': child infected during pregnancy', replacement = '', x = spec_output_tr$Source)) %>%
    dplyr::mutate(Source = gsub(pattern = ': child infected during breastfeeding', replacement = '', x = spec_output_tr$Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Started ART before', Source), 'Started ART before the pregnancy', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Started ART during', Source), 'Started ART during the pregnancy', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Started ART late', Source), 'Started ART late in the pregnancy', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Dropped off ART during pregnancy', Source), 'Dropped off ART', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Mother infected', Source), 'Mother infected', Source)) %>%
    dplyr::mutate(Source = ifelse(grepl('Did not receive ART', Source), 'Did not receive ART', Source)) %>%
    dplyr::mutate(spectrum  = spectrum)
  infections <- merge(frogger_output_infections, spec_output_infections, by = c('Source', 'trans_type', 'Year'))
  infections <- dplyr::mutate(infections, diff = frogger - spectrum)
  expect_true(all(abs(infections$diff) == 0 ))

})

test_that("Infections among children align", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)
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
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)

  spec_prev <- utils$offtrt

  hc1 <- dplyr::right_join((reshape2::melt(out$hc1_hiv_pop)), data.frame(Var1 = 1:7, cd4_cat = c("gte30", "26-30", "21-25", "16-20", "11-15", "5-10", "lte5")), by = "Var1")
  hc2 <- dplyr::right_join(data.frame(reshape2::melt(out$hc2_hiv_pop)), data.frame(Var1 = 1:6, cd4_cat = c("gte1000", "750-999", "500-749", "350-499", "200-349", "lte200")), by = "Var1")
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
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
  utils <- readRDS(test_path("testdata/child_test_utils.rds"))

  dp <- utils$dp
  pjnz <- utils$pjnz

  out <- run_model(parameters, "ChildModel", 1970:2030)

  spec_prev <- utils$ontrt
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
    dplyr::filter(year < 2030 & age < 15)

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("HIV related deaths among CLHIV not on ART align", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
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

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("HIV related deaths among CLHIV on ART align", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))
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

  expect_true(all(abs(dt$diff) < 5e-2))
})

test_that("Child model agrees when run through all years vs two parts vs single year runs", {
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))

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
  parameters <- read_parameters(test_path("testdata/child_parms.h5"))

  out_first_half_years <- run_model(parameters, "ChildModel", 1970:2000)
  expect_error(
    run_model_from_state(parameters, "ChildModel", get_time_slice(out_first_half_years, 31), 2000, 2000:2030),
    "Cannot output year '2000'. Output years must be later than simulation start year '2000'.",
    fixed = TRUE
  )
})
