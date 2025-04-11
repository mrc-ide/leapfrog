#' Run leapfrog model fit
#'
#' @param parameters Projection parameters
#' @param configuation The model configuration to run, see
#'   [list_model_configurations()] for available configurations
#' @param output_years Which years of the model to return from the simulation,
#'   defaults to all years from 1970 to 2030. Also used to control what years
#'   the simulation is run for. If output only 2030, simulation will be run
#'   from `projection_start_year` passed in the `parameters` list.
#'
#' @return List of model outputs
#' @export
run_model <- function(parameters,
                      configuration = "HivFullAgeStratification",
                      output_years = seq(1970, 2030)) {

  if (configuration == "HivCoarseAgeStratification") {
    parameters$hAG_SPAN <- parameters[["hAG_SPAN_coarse"]]
    parameters$cd4_initdist <- parameters[["cd4_initdist_coarse"]]
    parameters$cd4_prog <- parameters[["cd4_prog_coarse"]]
    parameters$cd4_mort <- parameters[["cd4_mort_coarse"]]
    parameters$art_mort <- parameters[["art_mort_coarse"]]
  } else if (configuration %in% c("HivFullAgeStratification", "ChildModel")) {
    parameters$hAG_SPAN <- parameters[["hAG_SPAN_full"]]
    parameters$cd4_initdist <- parameters[["cd4_initdist_full"]]
    parameters$cd4_prog <- parameters[["cd4_prog_full"]]
    parameters$cd4_mort <- parameters[["cd4_mort_full"]]
    parameters$art_mort <- parameters[["art_mort_full"]]
  }

  # convert indices to 0 based
  if ("artcd4elig_idx" %in% names(parameters)) {
    # integer type
    parameters[["artcd4elig_idx"]] <- parameters[["artcd4elig_idx"]] - 1L
  }
  if ("paed_art_elig_cd4" %in% names(parameters)) {
    # double type
    parameters[["paed_art_elig_cd4"]] <- parameters[["paed_art_elig_cd4"]] - 1
  }
  if ("t_ART_start" %in% names(parameters)) {
    parameters[["t_ART_start"]] <- parameters[["t_ART_start"]] - 1L
  }

  if (is_run_hiv_simulation(configuration)) {
    hTS <- dim(parameters[["art_mort"]])[[1]]
    parameters[["h_art_stage_dur"]] <- rep(0.5, hTS - 1)
  }

  parameters[["is_midyear_projection"]] <-
    parameters[["projection_period"]] == "midyear"

  run_base_model(parameters, configuration, output_years)
}

is_run_hiv_simulation <- function(configuration) {
  configuration %in% c("HivCoarseAgeStratification",
                       "HivFullAgeStratification",
                       "ChildModel")
}
