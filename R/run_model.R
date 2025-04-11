#' Run leapfrog model fit
#'
#' @param parameters Projection parameters
#' @param configuation The model configuration to run, see TODO
#' @param output_years Which years of the model to return from the simulation,
#'   defaults to all years from 1970 to 2030. Also used to control what years
#'   the simulation is run for. If output only 2030, simulation will be run
#'   from
#'
#' @return List of model outputs
#' @export
run_model <- function(parameters,
                      configuration = "HivFullAgeStratification",
                      output_years = seq(1970, 2030)) {

  assert_configuration_valid(configuration)

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
  configuration != "DemographicProjection"
}

assert_configuration_valid <- function(configuration) {
  valid_configs <- c("DemographicProjection",
                     "HivFullAgeStratification",
                     "HivCoarseAgeStratification",
                     "ChildModel")
  if (!(configuration %in% valid_configs)) {
    config_text <- paste(valid_configs, collapse = "', '")
    stop(sprintf(
      "Invalid configuration: '%s', must be one of '%s'.",
      configuration, config_text))
  }
}
