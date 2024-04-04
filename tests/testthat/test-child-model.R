test_that("child model can be run for all years", {
  input <- readRDS(test_path("testdata/child_parms.rds"))
  demp <- input$demp
  parameters <- input$proj
  parameters$ctx_effect <- 0.33
  parameters$laf <- 1
  parameters$paed_art_elig_age <- as.integer(parameters$paed_art_elig_age)
  parameters$mat_prev_input = rep(TRUE,61)

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

  ##Ensure that infections under one match
  dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }
  dp = input$dp
  u1_inf_spec <- dpsub("<NewInfantInfections MV>", 2, input$timedat.idx)
  ##this isn't true right now
  expect_true(abs(as.numeric(u1_inf_spec[7:length(input$timedat.idx)]) - colSums(out$p_infections[1,,7:length(input$timedat.idx)])) < 1e-1)



  ## All 10 as seeded with 100 infections which are distributed over 5 age
  ## groups and genders evenly
  expect_true(all(out$infections[1:5, , which(1970:2030 == 2000)] == 10))
  expect_true(all(out$hiv_population[1:5, , which(1970:2030 == 2000)] == 10))

  ## HIV population should be larger than zero in age 6 in 2001 because ageing
  ## is allowed
  expect_true(all(out$hc2_hiv_pop[1, 1, 0, , which(1970:2030 == 2001)] >= 0))

  ## HIV population should be larger than zero in CD4 categories after the
  ## highest as natural history now allowed
  expect_true(all(out$hc1_hiv_pop[2, 1, 3, , which(1970:2030 == 2001)] >= 0))


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
