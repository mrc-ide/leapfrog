test_that("Child model can be run for all years", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters

  expect_silent(out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE))
  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths", "p_hiv_pop",
      "p_hiv_pop_natural_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths", "hc1_hiv_pop", "hc2_hiv_pop",
      "hc1_art_pop", "hc2_art_pop",
      "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
      "hc1_art_aids_deaths", "hc2_art_aids_deaths", "hc_art_num", "hiv_births",
      "hc_art_total", "hc_art_init"
    )
  )

  ##Nothing should ever be negative
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
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)

  u1_inf_spec <- dpsub("<NewInfantInfections MV>", 2, input$timedat.idx)
  expect_true(all(abs(as.numeric(u1_inf_spec) - colSums(out$p_infections[1,,])) < 1e-3))
})

test_that("HIV related deaths among CLHIV not on ART align", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)
  tag.x ="<AIDSDeathsNoARTSingleAge MV>"
  start.id = 20898
  end.id = 21148
  timedat.idx = input$timedat.idx
  aids_deathsnoart <- array(as.numeric(unlist(dpsub(tag.x, 3:(end.id - start.id - 2), timedat.idx))), dim = c(length(3:(end.id - start.id - 2)),length(timedat.idx)))
  m = aids_deathsnoart[84:98,]
  f = aids_deathsnoart[166:180,]
  aids_deathsnoart <- array(0, dim = c(15,2,61))
  aids_deathsnoart[,1,] <- m
  aids_deathsnoart[,2,] <- f

  ##right now this is only working for the first year of ART, assuming its something with the timing on art
  hc1 <- apply(out$hc1_noart_aids_deaths, c(3:5), sum)
  hc2 <- apply(out$hc2_noart_aids_deaths, c(3:5), sum)
  hc <- array(0, dim = c(15,2,61))
  hc[1:5,,] <- hc1
  hc[6:15,,] <- hc2
  dt <- right_join(reshape2::melt(hc), reshape2::melt(aids_deathsnoart), by = c('Var1', 'Var2', 'Var3'))
  dt <- dt %>%
    mutate(age = Var1 - 1,
           sex = if_else(Var2 == 1, 'Male', 'Female'),
           year = Var3 + 1969,
           fr = value.x,
           spec = value.y) %>%
    select(age, sex, year, fr, spec)
  dt <- dt %>%
    mutate(diff = spec - fr)

  z = data.table(dt)
  expect_true(all(abs(dt$diff) < 1e-1))
})

test_that("CLHIV align", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)

  spec_prev <- input$offtrt
  hc1 = right_join(data.frame(reshape2::melt(out$hc1_hiv_pop)), data.frame(Var1 = 1:7, cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5')), by = 'Var1')
  hc2 = right_join(data.frame(reshape2::melt(out$hc2_hiv_pop)), data.frame(Var1 = 1:6, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200')), by = 'Var1')
  hc1 = right_join(hc1, data.frame(Var2 = 1:4, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+')), by = 'Var2')
  hc2 = right_join(hc2, data.frame(Var2 = 1:4, transmission = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+')), by = 'Var2')
  hc1 <- hc1 %>%
    mutate(age = Var3 - 1,
           sex = if_else(Var4 == 1, 'Male', 'Female'),
           year = Var5 + 1969,
           fr = value) %>%
    select(age, cd4_cat, transmission, sex, year, fr)
  hc2 <- hc2 %>%
    mutate(age = Var3 + 4,
           sex = if_else(Var4 == 1, 'Male', 'Female'),
           year = Var5 + 1969,
           fr = value) %>%
    select(age, cd4_cat,  transmission, sex, year, fr)
  hc <- rbind(hc1, hc2)
  hc <- hc %>%
    mutate(age = as.integer(age),
           cd4_cat = as.factor(cd4_cat),
           sex = as.factor(sex),
           year = as.integer(year),
           transmission = as.factor(transmission)) %>%
    group_by(age, cd4_cat, sex, year, transmission) %>%
    summarise(fr = sum(fr))
  dt <- right_join(hc, spec_prev, by = c('sex', 'age', 'cd4_cat', 'year', 'transmission'))
  dt <- dt %>%
    mutate(diff = pop - fr)
  x = data.table(dt)
  year.x = 2005;  x[year== year.x & age == 0 ] ; x[year== year.x & age == 1 & transmission == 'bf12+'];
  year.x = 2006;  x[year== year.x & age == 0 ] ; x[year== year.x & age == 1 & transmission == 'bf12+'];
  year.x = 2007;  x[year== year.x & age == 0 ] ; x[year== year.x & age == 1 & transmission == 'bf12+'];

  expect_true(all(abs(dt$diff) < 1e-1))
})

test_that("CLHIV on ART align", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)

  spec_prev <- input$ontrt
  spec_prev <- spec_prev %>%
    rename(time_art = transmission)
  hc1 = right_join(data.frame(reshape2::melt(out$hc1_art_pop)), data.frame(Var2 = 1:7, cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5')), by = ('Var2'))
  hc1 = right_join(hc1, data.frame(Var1 = 1:3, time_art = c('ARTlte5mo', 'ART6to12mo', 'ARTgte12mo')), by = ('Var1'))
  hc2 = right_join(data.frame(reshape2::melt(out$hc2_art_pop)), data.frame(Var2 = 1:6, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200')), by = 'Var2')
  hc2 = right_join(hc2, data.frame(Var1 = 1:3, time_art = c('ARTlte5mo', 'ART6to12mo', 'ARTgte12mo')), by = ('Var1'))
  hc1 <- hc1 %>%
    mutate(age = Var3 - 1,
           sex = if_else(Var4 == 1, 'Male', 'Female'),
           year = Var5 + 1969,
           fr = value) %>%
    select(age, cd4_cat, time_art, sex, year, fr)
  hc2 <- hc2 %>%
    mutate(age = Var3 + 4,
           sex = if_else(Var4 == 1, 'Male', 'Female'),
           year = Var5 + 1969,
           fr = value) %>%
    select(age, cd4_cat, time_art, sex, year, fr)
  hc <- rbind(hc1, hc2)
  hc <- hc %>%
    mutate(age = as.integer(age),
           cd4_cat = as.factor(cd4_cat),
           time_art = as.factor(time_art),
           sex = as.factor(sex),
           year = as.integer(year)) %>%
    group_by(age, cd4_cat, time_art, sex, year)
  dt <- right_join(hc, spec_prev, by = c('sex', 'age', 'cd4_cat', 'year', 'time_art'))
  dt <- dt %>%
    mutate(diff = pop - fr)

  expect_true(all(abs(dt$diff) < 1e-1))
})

#Not working
test_that("HIV related deaths among CLHIV on ART align", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)
  tag.x ="<AIDSDeathsARTSingleAge MV>"
  start.id = 20608
  end.id = 20858
  timedat.idx = input$timedat.idx
  aids_deathsart <- array(as.numeric(unlist(dpsub(tag.x, 3:(end.id - start.id - 2), timedat.idx))), dim = c(length(3:(end.id - start.id - 2)),length(timedat.idx)))
  m = aids_deathsart[84:98,]
  f = aids_deathsart[166:180,]
  aids_deathsart <- array(0, dim = c(15,2,61))
  aids_deathsart[,1,] <- m
  aids_deathsart[,2,] <- f

  ##right now this is only working for the first year of ART, assuming its something with the timing on art
  hc1 <- apply(out$hc1_art_aids_deaths, c(3:5), sum)
  hc2 <- apply(out$hc2_art_aids_deaths, c(3:5), sum)
  hc <- array(0, dim = c(15,2,61))
  hc[1:5,,] <- hc1
  hc[6:15,,] <- hc2
  dt <- right_join(reshape2::melt(hc), reshape2::melt(aids_deathsart), by = c('Var1', 'Var2', 'Var3'))
  dt <- dt %>%
    mutate(age = Var1 - 1,
           sex = if_else(Var2 == 1, 'Male', 'Female'),
           year = Var3 + 1969,
           fr = value.x,
           spec = value.y) %>%
    select(age, sex, year, fr, spec)
  dt <- dt %>%
    mutate(diff = spec - fr)

  out$hc1_art_aids_deaths[,,,,35]
  data.table(dt)[year == 2004 & age < 5]

  expect_true(all(abs(dt$diff) < 1))


})


