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
#' @return List of model outputs, where the last dimension of each element is
#'   time, e.g. `p_totpop` state variable has dimensions 81 x 2. If
#'   `output_years` specified has length 61 then the `p_totpop` output
#'   will have dimensions 81 x 2 x 61.
#' @examples
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "frogger", mustWork = TRUE)
#' parameters <- prepare_leapfrog_parameters(pjnz, use_coarse_age_groups = TRUE)
#' out <- run_model(parameters, "HivCoarseAgeStratification", 1970:2030)
#' @export
run_model <- function(parameters,
                      configuration = "HivFullAgeStratification",
                      output_years = seq(1970, 2030)) {
  parameters <- process_parameters_to_cpp(parameters)
  run_base_model(parameters, configuration, output_years)
}

#' Run leapfrog model fit from initial state
#'
#' @param parameters Projection parameters
#' @param configuration The model configuration to run, see
#'   [list_model_configurations()] for available configurations
#' @param initial_state The model will run from this initial state
#' @param simulation_start_year Start the model simulation from a particular year
#' @param output_years Which years of the model to return from the simulation,
#'   defaults to all years from 1970 to 2030. Also used to control what years
#'   the simulation is run for. If output only 2030, simulation will be run
#'   from `projection_start_year` passed in the `parameters` list.
#'
#' @return List of model outputs, where the last dimension of each element is
#'   time, e.g. `p_totpop` state variable has dimensions 81 x 2. If
#'   `output_years` specified has length 61 then the `p_totpop` output
#'   will have dimensions 81 x 2 x 61.
#' @examples
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "frogger", mustWork = TRUE)
#' parameters <- prepare_leapfrog_parameters(pjnz, use_coarse_age_groups = TRUE)
#' out_first_half_years <- run_model(parameters, "HivCoarseAgeStratification", 1970:2000)
#' out_second_half_years <- run_model_from_state(
#'   parameters,
#'   "HivCoarseAgeStratification",
#'   get_time_slice(out_first_half_years, 31),
#'   2000,
#'   2001:2030)
#' @export
run_model_from_state <- function(parameters,
                                 configuration,
                                 initial_state,
                                 simulation_start_year,
                                 output_years = seq(1970, 2030)) {
  parameters <- process_parameters_to_cpp(parameters)
  run_base_model_from_state(parameters, configuration, initial_state, simulation_start_year, output_years)
}

#' Run leapfrog model fit for a single year
#'
#' @param parameters Projection parameters
#' @param configuration The model configuration to run, see
#'   [list_model_configurations()] for available configurations
#' @param initial_state The model will run from this initial state
#' @param simulation_start_year Start the model simulation from this year
#'
#' @return List of model outputs without the last time dimension.
#'   This is different from [run_model_from_state()] and [run_model()]
#'   that do include the last time dimension. Since only the next time
#'   step is returned, dropping the time dimensions makes it easier to
#'   feed the returned list into the next single year model run. In
#'   contrast to [run_model_from_state()] the `p_totpop` output will
#'   have dimensions 81 x 2 not 81 x 2 x 61.
#' @examples
#' pjnz <- system.file(
#'   "pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
#'   package = "frogger", mustWork = TRUE)
#' parameters <- prepare_leapfrog_parameters(pjnz, use_coarse_age_groups = TRUE)
#' out_first_half_years <- run_model(parameters, "HivCoarseAgeStratification", 1970:2000)
#' prev_state <- get_time_slice(out_first_half_years, 31)
#' for (i in 2001:2029) {
#'   new_state <- run_model_single_year(parameters, "HivCoarseAgeStratification", prev_state, i)
#'   # Do things with new state, other processes, saving output etc.
#'   prev_state <- new_state
#' }
#' @export
run_model_single_year <- function(parameters,
                                  configuration,
                                  initial_state,
                                  simulation_start_year) {
  parameters <- process_parameters_to_cpp(parameters)
  run_base_model_single_year(parameters, configuration, initial_state, simulation_start_year)
}

#' Process parameters and convert from 1 based indexing in R to
#' 0 based indexing in C++. Also add in any defaults/extra parameters,
#' e.g. `h_art_stage_dur`.
#'
#' @param parameters The list of parameters to feed into the model.
#'
#' @return List of parameters with 0 based indexing and defaults.
#' @export
process_parameters_to_cpp <- function(parameters) {
  # convert indices to 0 based
  if ("artcd4elig_idx" %in% names(parameters)) {
    parameters[["artcd4elig_idx"]] <- parameters[["artcd4elig_idx"]] - 1L
  }
  if ("hc_art_elig_cd4" %in% names(parameters)) {
    parameters[["hc_art_elig_cd4"]] <- parameters[["hc_art_elig_cd4"]] - 1L
  }
  if ("t_ART_start" %in% names(parameters)) {
    parameters[["t_ART_start"]] <- parameters[["t_ART_start"]] - 1L
  }

  if ("artmx_timerr" %in% names(parameters)) {
    hTS <- dim(parameters[["artmx_timerr"]])[[1]]
    parameters[["h_art_stage_dur"]] <- rep(0.5, hTS - 1)
  }
  if (is.null(parameters[["hts_per_year"]])) {
    parameters[["hts_per_year"]] <- 10L
  }

  parameters
}

#' Process parameters and convert from 0 based indexing for C++ to
#' 1 based indexing in R.
#'
#' @param parameters List of parameters.
#'
#' @return List of parameters with 1 based indexing.
#' @export
process_parameters_to_r <- function(parameters) {
  # convert indices to 0 based
  if ("artcd4elig_idx" %in% names(parameters)) {
    parameters[["artcd4elig_idx"]] <- parameters[["artcd4elig_idx"]] + 1L
  }
  if ("hc_art_elig_cd4" %in% names(parameters)) {
    parameters[["hc_art_elig_cd4"]] <- parameters[["hc_art_elig_cd4"]] + 1L
  }
  if ("t_ART_start" %in% names(parameters)) {
    parameters[["t_ART_start"]] <- parameters[["t_ART_start"]] + 1L
  }

  parameters
}

#' Slice a single year from model state
#'
#' @param state The model state with time dimension
#' @param index The index of the time step you want to extract
#'
#' @return List of model outputs for the specified time step. Can be used as
#'   input state for [run_model_from_state()] and [run_model_single_year()].
#'   All outputs will have 1 fewer dimension than input state.
#' @export
get_time_slice <- function(state, index) {
  last_ind <- function(x) {
    nd <- length(dim(x))
    inds <- rep(alist(,)[1], nd)
    inds[nd] <- index
    ret <- do.call(`[`, c(list(x), inds))

    # R drops dimension by default when extracting a 1D
    # vector which causes waldo compare to fail so we
    # manually compute the dimension for those fields.
    if (is.null(dim(ret)) && length(ret) > 1) {
      dim(ret) <- length(ret)
    }
    ret
  }

  lapply(state, last_ind)
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
