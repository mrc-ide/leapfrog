import numpy as np
from leapfrog import (  # type: ignore[attr-defined]
    Art,
    BaseModelChildParameters,
    BaseModelParameters,
    BaseModelState,
    ChildModelState,
    Demography,
    Incidence,
    NaturalHistory,
    Options,
    Parameters,
    State,
    project_single_year_cpp,
    run_model_cpp,
    set_initial_state_cpp,
)

from leapfrog_py.utils import assert_fortran_order, assert_keys


def set_initial_state(
    parameters: dict[str, np.ndarray], state: dict[str, np.ndarray]
):
    assert_fortran_order(state)
    params = _initialise_params(parameters, 10)
    initial_state = _initialise_state(0, state, params)
    set_initial_state_cpp(initial_state, params)


def project_single_year(
    i: int,
    parameters: dict[str, np.ndarray],
    state: dict[str, np.ndarray],
    *,
    hts_per_year: int = 10,
):
    assert_fortran_order(state)
    if i < 1:
        msg = "Year must be greater than 0, 0th year is the initial state. See 'set_initial_state'."
        raise ValueError(msg)

    params = _initialise_params(parameters, hts_per_year)
    state_prev = _initialise_state(i - 1, state, params)
    state_now = _initialise_state(i, state, params)

    project_single_year_cpp(i, params, state_prev, state_now)


def run_leapfrog(
    parameters: dict[str, np.ndarray], *, sim_years=None, hts_per_year=10
) -> dict[str, np.ndarray]:
    if sim_years is None:
        sim_years = np.arange(1970, 2031)

    time_steps = len(sim_years)
    save_steps = range(time_steps)
    params = _initialise_params(parameters, hts_per_year)

    return run_model_cpp(params, time_steps, save_steps)


def _initialise_params(
    params: dict[str, np.ndarray], hts_per_year: int
) -> Parameters:
    required_input = [
        "basepop",
        "survival_probability",
        "net_migration",
        "age_specific_fertility_rate",
        "births_sex_prop",
        "incidence_rate",
        "incidence_age_rate_ratio",
        "incidence_sex_rate_ratio",
        "cd4_mortality_full",
        "cd4_progression_full",
        "cd4_initial_distribution_full",
        "idx_hm_elig",
        "art_mortality_rate_full",
        "art_mortality_time_rate_ratio",
        "art_dropout",
        "adults_on_art",
        "adults_on_art_is_percent",
        "h_art_stage_dur",
    ]
    assert_keys(params, required_input)
    assert_fortran_order(params)

    demography = Demography(
        params["basepop"],
        params["survival_probability"],
        params["net_migration"],
        params["age_specific_fertility_rate"],
        params["births_sex_prop"],
    )
    incidence = Incidence(
        params["incidence_rate"],
        params["incidence_age_rate_ratio"],
        params["incidence_sex_rate_ratio"],
    )
    nat_history = NaturalHistory(
        params["cd4_mortality_full"],
        params["cd4_progression_full"],
        params["cd4_initial_distribution_full"],
        1,
    )
    art = Art(
        params["idx_hm_elig"],
        params["art_mortality_rate_full"],
        params["art_mortality_time_rate_ratio"],
        params["art_dropout"],
        params["adults_on_art"],
        params["adults_on_art_is_percent"],
        params["h_art_stage_dur"],
        1,
    )
    # TODO: Move this into C++, 30 is time ART start and 66 is no of HIV age groups which should come from the
    # state space in C++
    options = Options(hts_per_year, 30, 66)
    base_model_params = BaseModelParameters(
        options, demography, incidence, nat_history, art
    )
    child_model_params = BaseModelChildParameters()
    return Parameters(base_model_params, child_model_params)


def _initialise_state(
    time_step: int, state: dict[str, np.ndarray], params
) -> State:
    year_state = {k: np.atleast_1d(v[..., time_step]) for k, v in state.items()}
    base_model_state = BaseModelState(
        p_total_pop=year_state["p_total_pop"],
        births=year_state["births"],
        p_total_pop_natural_deaths=year_state["p_total_pop_natural_deaths"],
        p_hiv_pop=year_state["p_hiv_pop"],
        p_hiv_pop_natural_deaths=year_state["p_hiv_pop_natural_deaths"],
        h_hiv_adult=year_state["h_hiv_adult"],
        h_art_adult=year_state["h_art_adult"],
        h_hiv_deaths_no_art=year_state["h_hiv_deaths_no_art"],
        p_infections=year_state["p_infections"],
        h_hiv_deaths_art=year_state["h_hiv_deaths_art"],
        h_art_initiation=year_state["h_art_initiation"],
        p_hiv_deaths=year_state["p_hiv_deaths"],
    )
    # TODO: Initialise child model state from the data, and remove params from the function signature
    return State(base_model_state, ChildModelState(params))
