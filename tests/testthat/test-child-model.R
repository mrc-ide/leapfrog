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
      "hc1_art_aids_deaths", "hc2_art_aids_deaths", "hc_art_num", "hiv_births"
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


test_that("HIV related deaths among children align", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)

  start.id = 21292
  end.id = 21458
  aidsdeaths <- array(as.numeric(unlist(dpsub("<AidsDeathsByAge MV2>"  , 3:(end.id - start.id - 2), timedat.idx))), dim = c(length(3:(end.id - start.id - 2)),length(timedat.idx)))
  m = aidsdeaths[1:15,]
  f = aidsdeaths[82:96,]
  aidsdeaths <- array(0, dim = c(15,2,61))
  aidsdeaths[,1,] <- m
  aidsdeaths[,2,] <- f
  aidsdeaths[,,which(1970:2030 == 2002)]-
  out$p_hiv_deaths[1:15,,which(1970:2030 == 2002)]
  expect_true(all(abs(aidsdeaths - out$p_hiv_deaths[1:15,,]) < 1e-1))
})

test_that("CLHIV align", {
  input <- setup_childmodel(testinput = "testdata/child_parms.rds")
  demp = input$demp
  parameters = input$parameters
  dp = input$dp
  pjnz = input$pjnz

  out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE)

  spec_prev <- input$pop1
  hc1 = right_join(data.frame(reshape2::melt(out$hc1_hiv_pop)), data.frame(Var1 = 1:7, cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5')), by = 'Var1')
  hc2 = right_join(data.frame(reshape2::melt(out$hc2_hiv_pop)), data.frame(Var1 = 1:6, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200')), by = 'Var1')
  hc1 <- hc1 %>%
    mutate(age = Var3 - 1,
           sex = if_else(Var4 == 1, 'Male', 'Female'),
           year = Var5 + 1969,
           fr = value) %>%
    select(age, cd4_cat, sex, year, fr)
  hc2 <- hc2 %>%
    mutate(age = Var3 + 4,
           sex = if_else(Var4 == 1, 'Male', 'Female'),
           year = Var5 + 1969,
           fr = value) %>%
    select(age, cd4_cat, sex, year, fr)
  hc <- rbind(hc1, hc2)
  hc <- hc %>%
    mutate(age = as.integer(age),
           cd4_cat = as.factor(cd4_cat),
           sex = as.factor(sex),
           year = as.integer(year)) %>%
    group_by(age, cd4_cat, sex, year) %>%
    summarise(fr = sum(fr))
  dt <- right_join(hc, spec_prev, by = c('sex', 'age', 'cd4_cat', 'year'))
  dt <- dt %>%
    mutate(diff = pop - fr)
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

