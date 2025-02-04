setup_PaediatricModel <- function(testinput = "testdata/child_parms.rds") {
  input <- readRDS(testthat::test_path(testinput))
  demp <- input$demp
  parameters <- input$proj

  #### TO DO: move this into the frogger version of spectrum-inputs.R
  parameters$ctx_effect <- 0.33
  parameters$laf <- 1
  parameters$paed_art_elig_age <- as.integer(parameters$paed_art_elig_age)
  parameters$mat_prev_input <- rep(TRUE, 61)
  pmtct_new <- array(0, dim = c(7, 61), dimnames = list(pmtct = c("Option A", "Option B", "SDNVP", "Dual ARV", "Option B+: before pregnancy", "Option B+: >4 weeks", "Option B+: <4 weeks")))
  ## pick out which ones were inserted as numbers
  pmtct_new[, which(colSums(parameters$pmtct)[, 1] > 0)] <- parameters$pmtct[, (which(colSums(parameters$pmtct)[, 1] > 0)), 1]
  ## pick out which ones were inserted as percent
  pmtct_new[, which(colSums(parameters$pmtct)[, 1] == 0)] <- parameters$pmtct[, which(colSums(parameters$pmtct)[, 1] == 0), 2]
  parameters$pmtct <- pmtct_new

  return(list(
    dp = input$dp, demp = demp, parameters = parameters,
    pjnz = input$pjnz, timedat.idx = input$timedat.idx,
    pop1 = input$pop1_outputs,
    ontrt = input$on_treatment,
    offtrt = input$off_trt,
    deaths_noart = input$deaths_noart,
    deaths_art = input$deaths_art
  ))
}

dpsub <- function(tag, rows, cols, tagcol = 1) {
  dp[which(dp[, tagcol] == tag) + rows, cols]
}
