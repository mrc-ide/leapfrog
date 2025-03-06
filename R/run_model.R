#' Run leapfrog model fit
#'
#' @param data Input data
#' @param parameters Projection parameters
#' @param sim_years Simulation years to run model for default 1970:2030
#' @param hts_per_year Number of HIV time steps per year, default 10
#' @param output_steps Which sim years to output for default same as sim_years
#' @param run_hiv_simulation If TRUE then runs HIV simulation
#' @param hiv_age_stratification The age stratification for HIV population,
#'  "coarse" or "full"
#' @param run_child_model If TRUE then run the child model
#'
#' @return List of model outputs
#' @export
run_model <- function(data, parameters, sim_years,
                      hts_per_year,
                      output_steps = NULL,
                      run_hiv_simulation = TRUE,
                      hiv_age_stratification = "full",
                      run_child_model = TRUE) {
  if (is.null(sim_years)) {
    sim_years <- 1970:2030
  }
  if (is.null(output_steps)) {
    output_steps <- seq_along(sim_years)
  } else {
    ## We want users to think in terms of years, so interface has years
    ## But for running the C++ loop we want to report out based on what
    ## iteration step we're at so convert from years to index
    invalid_steps <- output_steps[!(output_steps %in% sim_years)]
    if (any(invalid_steps)) {
      out_str <- paste(paste0("'", invalid_steps, "'"), collapse = ", ")
      stop(sprintf(
        "Invalid output %s %s. Can only output one of the simulation years.",
        ngettext(length(invalid_steps), "step", "steps"), out_str
      ))
    }
    output_steps <- which(sim_years %in% output_steps)
  }

  assert_enum(hiv_age_stratification, c("full", "coarse", "none"))
  if (!run_hiv_simulation) {
    model_variant <- "DemographicProjection"
  } else if (hiv_age_stratification == "full" && !run_child_model) {
    model_variant <- "HivFullAgeStratification"
  } else if (hiv_age_stratification == "coarse" && !run_child_model) {
    model_variant <- "HivCoarseAgeStratification"
  } else if (hiv_age_stratification == "full" && run_child_model) {
    model_variant <- "ChildModel"
  } else if (hiv_age_stratification == "coarse" && run_child_model) {
    stop("Cannot run child model with coarse age stratification")
  } else if (hiv_age_stratification == "none" && run_hiv_simulation) {
    stop("Cannot run HIV simulation with unspecified age stratification")
  }
  if (hiv_age_stratification == "full") {
    parameters$hAG_SPAN <- parameters[["hAG_SPAN_full"]]
    parameters$cd4_initdist <- parameters[["cd4_initdist_full"]]
    parameters$cd4_prog <- parameters[["cd4_prog_full"]]
    parameters$cd4_mort <- parameters[["cd4_mort_full"]]
    parameters$art_mort <- parameters[["art_mort_full"]]
  } else if (hiv_age_stratification == "coarse") {
    parameters$hAG_SPAN <- parameters[["hAG_SPAN_coarse"]]
    parameters$cd4_initdist <- parameters[["cd4_initdist_coarse"]]
    parameters$cd4_prog <- parameters[["cd4_prog_coarse"]]
    parameters$cd4_mort <- parameters[["cd4_mort_coarse"]]
    parameters$art_mort <- parameters[["art_mort_coarse"]]
  }
  data <- c(data, parameters)

  if (is.null(hts_per_year)) {
    hts_per_year <- 10
  }

  # convert indices to 0 based
  if ("artcd4elig_idx" %in% names(data)) {
    # integer type
    data[["artcd4elig_idx"]] <- data[["artcd4elig_idx"]] - 1L
  }
  if ("paed_art_elig_cd4" %in% names(data)) {
    # double type
    data[["paed_art_elig_cd4"]] <- data[["paed_art_elig_cd4"]] - 1
  }

  if (run_hiv_simulation) {
    hTS <- dim(data[["art_mort"]])[[1]]
    data[["h_art_stage_dur"]] <- rep(0.5, hTS - 1)
  }

  run_base_model(data, model_variant, length(sim_years), hts_per_year, output_steps - 1, data[["projection_period"]] == "midyear", data[["t_ART_start"]] - 1)
}
