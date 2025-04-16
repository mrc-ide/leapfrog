#' Run leapfrog model fit
#'
#' @param parameters Projection parameters
#' @param configuration The model configuration to run, see
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
  parameters <- process_parameters(parameters, configuration)

  run_base_model(parameters, configuration, output_years)
}

#' Run leapfrog model fit from initial state
#'
#' @param parameters Projection parameters
#' @param configuration The model configuration to run, see
#'   [list_model_configurations()] for available configurations
#' @param initial_state The model will run from this initial state
#' @param start_from_year Start the model simulation from a particular year
#' @param output_years Which years of the model to return from the simulation,
#'   defaults to all years from 1970 to 2030. Also used to control what years
#'   the simulation is run for. If output only 2030, simulation will be run
#'   from `projection_start_year` passed in the `parameters` list.
#'
#' @return List of model outputs
#' @export
run_model_from_state <- function(parameters,
                                 configuration,
                                 initial_state,
                                 start_from_year,
                                 output_years = seq(1970, 2030)) {
  parameters <- process_parameters(parameters, configuration)

  run_base_model_from_state(parameters, configuration, initial_state, start_from_year, output_years)
}


process_parameters <- function(parameters, configuration) {
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
  if (is.null(parameters[["hts_per_year"]])) {
    parameters[["hts_per_year"]] <- 10L
  }

  parameters
}

is_run_hiv_simulation <- function(configuration) {
  configuration %in% c("HivCoarseAgeStratification",
                       "HivFullAgeStratification",
                       "ChildModel")
}

get_last_time_slice <- function(dat) {
  dims <- dim(dat[[1]])
  last_t_index <- dims[[length(dims)]]

  last_ind <- function(x) {
    nd <- length(dim(x))
    inds <- rep(alist(,)[1], nd)
    inds[nd] <- last_t_index
    do.call(`[`, c(list(x), inds))
  }

  lapply(dat, last_ind)
}

concat_on_time_dim <- function(dat1, dat2) {
  x <- lapply(names(dat1), function(name) {
    con <- abind::abind(dat1[[name]], dat2[[name]])
    dimnames(con) <- NULL
    con
  })
  names(x) <- names(dat1)
  x
}
