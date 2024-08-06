#' Run leapfrog model fit
#'
#' @param data Input data
#' @param parameters Projection parameters
#' @param sim_years Simulation years to run model for default 1970:2030
#' @param hts_per_year Number of HIV time steps per year, default 10
#' @param output_steps Which sim years to output for default same as sim_years
#' @param hiv_age_stratification The age stratification for HIV population,
#'  "coarse" or "full"
#' @param run_child_model If TRUE then run the child model
#'
#' @return List of model outputs
#' @export
run_model <- function(data, parameters, sim_years,
                      hts_per_year,
                      output_steps = NULL,
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

  assert_enum(hiv_age_stratification, c("full", "coarse"))
  if (hiv_age_stratification == "full" && !run_child_model) {
    model_variant <- "BaseModelFullAgeStratification"
  } else if (hiv_age_stratification == "coarse" && !run_child_model) {
    model_variant <- "BaseModelCoarseAgeStratification"
  } else if (hiv_age_stratification == "full" && run_child_model) {
    model_variant <- "ChildModel"
  } else if (hiv_age_stratification == "coarse" && run_child_model) {
    stop("Cannot run child model with coarse age stratification")
  }
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


  proj_years <- transform_simulation_years(data, sim_years)
  save_steps <- transform_output_steps(output_steps)
  hiv_steps <- transform_hts_per_year(hts_per_year)

  run_base_model(data, proj_years, hiv_steps, save_steps, model_variant)
}

transform_simulation_years <- function(data, sim_years) {
  sx <- data$Sx
  max_sim_years <- dim(sx)[3]
  if (is.null(sim_years)) {
    return(max_sim_years)
  }
  n_years <- length(sim_years)
  if (n_years > max_sim_years) {
    stop(paste0("No of years > max years of ", max_sim_years))
  }
  n_years
}

transform_output_steps <- function(output_steps) {
  output_steps - 1
}

transform_hts_per_year <- function(hts_per_year) {
  if (is.null(hts_per_year)) {
    return(10)
  }
  hts_per_year
}

test_tmb <- function(data, parameters, sim_years,
                     hts_per_year,
                     output_steps = NULL,
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

  assert_enum(hiv_age_stratification, c("full", "coarse"))
  if (hiv_age_stratification == "full" && !run_child_model) {
    model_variant <- "BaseModelFullAgeStratification"
  } else if (hiv_age_stratification == "coarse" && !run_child_model) {
    model_variant <- "BaseModelCoarseAgeStratification"
  } else if (hiv_age_stratification == "full" && run_child_model) {
    model_variant <- "ChildModel"
  } else if (hiv_age_stratification == "coarse" && run_child_model) {
    stop("Cannot run child model with coarse age stratification")
  }
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


  proj_years <- transform_simulation_years(data, sim_years)
  save_steps <- transform_output_steps(output_steps)
  hiv_steps <- transform_hts_per_year(hts_per_year)

  # TMB specific stuff
  names_data <- names(data)
  for (name in names_data) {
    if (length(data[[name]]) > 1 && !is.array(data[[name]])) {
      data[[name]] <- array(data[[name]])
    }
  }

  parameter_tmb <- list(incidinput_scalar = 0.5)

  data_tmb <- list(
    data_vars = data,
    t_ART_start = data$t_ART_start,
    model_variant = model_variant,
    proj_years = proj_years,
    hiv_steps = hiv_steps,
    save_steps = save_steps,
    basepop = data$basepop,
    sd = 10000,
    incidinput_data = data$incidinput
  )

  obj <- TMB::MakeADFun(data_tmb, parameter_tmb, DLL = "frogger_TMB")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  list(
    obj = obj,
    opt = opt
  )
}
