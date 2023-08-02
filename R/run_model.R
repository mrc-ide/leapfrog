#' Run leapfrog model fit
#'
#' @param data Input data
#' @param parameters Projection parameters
#' @param sim_years Simulation years to run model for
#' @param hts_per_year Number of HIV time steps per year
#' @param output_steps Which sim years to output for
#' @param hiv_age_stratification The age stratification for HIV population,
#'  "coarse" or "full"
#' @param run_child_model If TRUE then run the child model
#'
#' @return List of model outputs
#' @export
run_model <- function(data, parameters, sim_years,
                      hts_per_year, output_steps,
                      hiv_age_stratification = "full",
                      run_child_model = TRUE) {
  assert_enum(hiv_age_stratification, c("full", "coarse"))
  if (hiv_age_stratification == "full") {
    parameters$hAG_SPAN <- parameters[["hAG_SPAN_full"]]
    parameters$cd4_initdist <- parameters[["cd4_initdist_full"]]
    parameters$cd4_prog <- parameters[["cd4_prog_full"]]
    parameters$cd4_mort <- parameters[["cd4_mort_full"]]
    parameters$art_mort <- parameters[["art_mort_full"]]
  } else {
    parameters$hAG_SPAN <- parameters[["hAG_SPAN_coarse"]]
    parameters$cd4_initdist <- parameters[["cd4_initdist_coarse"]]
    parameters$cd4_prog <- parameters[["cd4_prog_coarse"]]
    parameters$cd4_mort <- parameters[["cd4_mort_coarse"]]
    parameters$art_mort <- parameters[["art_mort_coarse"]]
  }
  data <- c(data, parameters)
  run_base_model(data, sim_years, hts_per_year, output_steps,
                 hiv_age_stratification, run_child_model)
}
