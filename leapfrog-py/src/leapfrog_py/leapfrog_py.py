import numpy as np
from leapfrog import (  # type: ignore[attr-defined]
    BaseModelState,
    ChildModelState,
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
    _validate_params(parameters)
    initial_state = _initialise_state(0, state)
    set_initial_state_cpp(initial_state, parameters)


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

    state_prev = _initialise_state(i - 1, state)
    state_now = _initialise_state(i, state)

    project_single_year_cpp(i, hts_per_year, parameters, state_prev, state_now)


def run_leapfrog(
    parameters: dict[str, np.ndarray], *, sim_years=None, hts_per_year=10
) -> dict[str, np.ndarray]:
    _validate_params(parameters)
    if sim_years is None:
        sim_years = np.arange(1970, 2031)

    time_steps = len(sim_years)
    save_steps = range(time_steps)

    return run_model_cpp(parameters, time_steps, hts_per_year, save_steps)


def _validate_params(params: dict[str, np.ndarray]) -> None:
    required_input = [
        "basepop",
        "Sx",
        "netmigr_adj",
        "asfr",
        "births_sex_prop",
        "incidinput",
        "incrr_age",
        "incrr_sex",
        "cd4_mort",
        "cd4_prog",
        "cd4_initdist",
        "artcd4elig_idx",
        "art_mort",
        "artmx_timerr",
        "art_dropout",
        "art15plus_num",
        "art15plus_isperc",
        "h_art_stage_dur",
        "t_ART_start",
        "scale_cd4_mort",
        "art_alloc_mxweight"
    ]
    assert_keys(params, required_input)
    assert_fortran_order(params)


def _initialise_state(time_step: int, state: dict[str, np.ndarray]) -> State:
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
        p_hiv_deaths=year_state["p_hiv_deaths"]
    )
    child_model_state = ChildModelState()
    return State(base_model_state, child_model_state)


def _write_state(
    time_step: int, state: State, full_state: dict[str, np.ndarray]
):
    # TODO: Code generate this, possibly on the C++ side. But hopefully it will be obsolete if we can
    # pass in a numpy view of the state by reference, to make this copying unnecessary.
    full_state["p_total_pop"][..., time_step] = state.base.p_total_pop
    full_state["births"][..., time_step] = state.base.births
    full_state["p_total_pop_natural_deaths"][
        ..., time_step
    ] = state.base.p_total_pop_natural_deaths
    full_state["p_hiv_pop"][..., time_step] = state.base.p_hiv_pop
    full_state["p_hiv_pop_natural_deaths"][
        ..., time_step
    ] = state.base.p_hiv_pop_natural_deaths
    full_state["h_hiv_adult"][..., time_step] = state.base.h_hiv_adult
    full_state["h_art_adult"][..., time_step] = state.base.h_art_adult
    full_state["h_hiv_deaths_no_art"][
        ..., time_step
    ] = state.base.h_hiv_deaths_no_art
    full_state["p_infections"][..., time_step] = state.base.p_infections
    full_state["h_hiv_deaths_art"][..., time_step] = state.base.h_hiv_deaths_art
    full_state["h_art_initiation"][..., time_step] = state.base.h_art_initiation
    full_state["p_hiv_deaths"][..., time_step] = state.base.p_hiv_deaths
