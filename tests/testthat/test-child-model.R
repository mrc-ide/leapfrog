test_that("child model can be run for all years", {
  input <- readRDS(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$proj
  parameters$ctx_effect <- 0.33
  parameters$laf <- 1
  parameters$paed_art_elig_age <- as.integer(parameters$paed_art_elig_age)
  parameters$mat_prev_input = rep(TRUE,61)
  parameters$pmtct <- parameters$pmtct[,,2]

  expect_silent(out <- run_model(demp, parameters, NULL, NULL, 0:60, run_child_model = TRUE))
  expect_setequal(
    names(out),
    c(
      "p_total_pop", "births", "p_total_pop_natural_deaths", "p_hiv_pop",
      "p_hiv_pop_natural_deaths", "h_hiv_adult", "h_art_adult",
      "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
      "h_art_initiation", "p_hiv_deaths", "hc1_hiv_pop", "hc2_hiv_pop",
      "hc1_art_pop", "hc2_art_pop", "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
      "hc1_art_aids_deaths", "hc2_art_aids_deaths", "hc_art_num", "hiv_births"
    )
  )

  expect_true(all(out$hc1_hiv_pop >= 0))

  ##load in DP functions
  pjnz = input$pjnz
  dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }
  dp = input$dp

  ##Ensure that infections under one match
  u1_inf_spec <- dpsub("<NewInfantInfections MV>", 2, input$timedat.idx)
  expect_true(all(abs(as.numeric(u1_inf_spec[7:length(input$timedat.idx)]) - colSums(out$p_infections[1,,7:length(input$timedat.idx)])) < 1e1))

  #Ensure that deaths align
  aidsdeaths <- array(as.numeric(unlist(dpsub("<AidsDeathsByAge MV2>"  , 3:(end.id - start.id - 2), timedat.idx))), dim = c(length(3:(end.id - start.id - 2)),length(timedat.idx)))
  m = aidsdeaths[1:15,]
  f = aidsdeaths[82:96,]
  aidsdeaths <- array(0, dim = c(15,2,61))
  aidsdeaths[,1,] <- m
  aidsdeaths[,2,] <- f
  expect_true(all(abs(aidsdeaths - out$p_hiv_deaths[1:15,,]) < 1e-1))

  #Ensure that prevalence aligns
  spec_prev <- input$pop1_outputs
  hc1 = merge(data.table(melt(out$hc1_hiv_pop)), data.table(Var1 = 1:7, cd4_cat = c('gte30', '26-30', '21-25', '16-20', '11-5', '5-10', 'lte5')), by = 'Var1')
  hc2 = merge(data.table(melt(out$hc2_hiv_pop)), data.table(Var1 = 1:6, cd4_cat = c('gte1000', '750-999', '500-749', '350-499', '200-349','lte200')), by = 'Var1')
  hc1 <- hc1[Var2 == 1,.(age = Var3 - 1, cd4_cat, sex = ifelse(Var4 == 1, 'Male', 'Female'), year = Var5 + 1969, fr = value)]
  hc2 <- hc2[Var2 == 1,.(age = Var3 + 4, cd4_cat, sex = ifelse(Var4 == 1, 'Male', 'Female'), year = Var5 + 1969, fr = value)]
  hc <- rbind(hc1, hc2)
  dt <- merge(hc, spec_prev, by = c('sex', 'age', 'cd4_cat', 'year'))
  dt[,diff := pop - fr]
  expect_true(all(abs(dt$diff) < 1e-1))


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
