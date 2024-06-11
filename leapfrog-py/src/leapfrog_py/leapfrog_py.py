import numpy as np
from leapfrog import (
    Parameters,
    Demography,
    Incidence,
    NaturalHistory,
    Art,
    Options,
    BaseModelParameters,
    BaseModelChildParameters,
    State,
    BaseModelState,
    ChildModelState,
    run_model_cpp,
    project_single_year_cpp,
    set_initial_state_cpp
)
from leapfrog_py.utils import assert_keys, assert_fortran_order


def set_initial_state(parameters: dict[str, np.ndarray], state: dict[str, np.ndarray]):
    assert_fortran_order(state)
    params = _initialise_params(parameters, 10)
    initial_state = _initialise_state(0, state, params)
    set_initial_state_cpp(initial_state, params)
    _write_state(0, initial_state, state)


def project_single_year(i: int, parameters: dict[str, np.ndarray], state: dict[str, np.ndarray], *,
                        hts_per_year: int = 10):
    assert_fortran_order(state)
    if i < 1:
        raise ValueError("Year must be greater than 0, 0th year is the initial state. See 'set_initial_state'.")

    params = _initialise_params(parameters, hts_per_year)
    state_prev = _initialise_state(i - 1, state, params)
    state_now = _initialise_state(i, state, params)

    project_single_year_cpp(i, params, state_prev, state_now)

    _write_state(i, state_now, state)


def run_leapfrog(parameters: dict[str, np.ndarray], *, sim_years=None, hts_per_year=10) -> dict[str, np.ndarray]:
    if sim_years is None:
        sim_years = np.arange(1970, 2031)

    time_steps = len(sim_years)
    save_steps = range(time_steps)
    params = _initialise_params(parameters, hts_per_year)

    return run_model_cpp(params, time_steps, save_steps)


def _initialise_params(params: dict[str, np.ndarray], hts_per_year: int) -> Parameters:
    required_input = ['basepop', 'survival_probability', 'net_migration', 'age_specific_fertility_rate',
                      'births_sex_prop', 'incidence_rate', 'relative_risk_age', 'relative_risk_sex',
                      'cd4_mortality_full', 'cd4_progression_full', 'cd4_initial_distribution_full', 'idx_hm_elig',
                      'art_mortality_rate_full', 'art_mortality_time_rate_ratio', 'art_dropout', 'adult_on_art',
                      'adults_on_art_is_percent', 'h_art_stage_dur']
    assert_keys(params, required_input)
    assert_fortran_order(params)

    demography = Demography(
        params['basepop'],
        params['survival_probability'],
        params['net_migration'],
        params['age_specific_fertility_rate'],
        params['births_sex_prop']
    )
    incidence = Incidence(
        params['incidence_rate'],
        params['relative_risk_age'],
        params['relative_risk_sex']
    )
    nat_history = NaturalHistory(
        params['cd4_mortality_full'],
        params['cd4_progression_full'],
        params['cd4_initial_distribution_full'],
        1
    )
    art = Art(
        params['idx_hm_elig'],
        params['art_mortality_rate_full'],
        params['art_mortality_time_rate_ratio'],
        params['art_dropout'],
        params['adult_on_art'],
        params['adults_on_art_is_percent'],
        params['h_art_stage_dur'],
        1
    )
    options = Options(hts_per_year, 30, 66)
    base_model_params = BaseModelParameters(options, demography, incidence, nat_history, art)
    child_model_params = BaseModelChildParameters()
    return Parameters(base_model_params, child_model_params)


def _initialise_state(time_step: int, state: dict[str, np.ndarray], params) -> State:
    state = {k: v[..., time_step] for k, v in state.items()}
    base_model_state = BaseModelState(
        p_total_pop=state['p_total_pop'],
        births=state['births'].item(),
        p_total_pop_natural_deaths=state['p_total_pop_natural_deaths'],
        p_hiv_pop=state['p_hiv_pop'],
        p_hiv_pop_natural_deaths=state['p_hiv_pop_natural_deaths'],
        h_hiv_adult=state['h_hiv_adult'],
        h_art_adult=state['h_art_adult'],
        h_hiv_deaths_no_art=state['h_hiv_deaths_no_art'],
        p_infections=state['p_infections'],
        h_hiv_deaths_art=state['h_hiv_deaths_art'],
        h_art_initiation=state['h_art_initiation'],
        p_hiv_deaths=state['p_hiv_deaths']
    )
    # TODO: Initialise child model state from the data, and remove params from the function signature
    return State(base_model_state, ChildModelState(params))


def _write_state(time_step: int, state: State, full_state: dict[str, np.ndarray]):
    # TODO: Code generate this, possibly on the C++ side. But hopefully it will be obsolete if we can
    # pass in a numpy view of the state by reference, to make this copying unnecessary.
    full_state['p_total_pop'][..., time_step] = state.base.p_total_pop
    full_state['births'][..., time_step] = state.base.births
    full_state['p_total_pop_natural_deaths'][..., time_step] = state.base.p_total_pop_natural_deaths
    full_state['p_hiv_pop'][..., time_step] = state.base.p_hiv_pop
    full_state['p_hiv_pop_natural_deaths'][..., time_step] = state.base.p_hiv_pop_natural_deaths
    full_state['h_hiv_adult'][..., time_step] = state.base.h_hiv_adult
    full_state['h_art_adult'][..., time_step] = state.base.h_art_adult
    full_state['h_hiv_deaths_no_art'][..., time_step] = state.base.h_hiv_deaths_no_art
    full_state['p_infections'][..., time_step] = state.base.p_infections
    full_state['h_hiv_deaths_art'][..., time_step] = state.base.h_hiv_deaths_art
    full_state['h_art_initiation'][..., time_step] = state.base.h_art_initiation
    full_state['p_hiv_deaths'][..., time_step] = state.base.p_hiv_deaths
