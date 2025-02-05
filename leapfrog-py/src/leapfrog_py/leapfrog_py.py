import numpy as np
from leapfrog import (  # type: ignore[attr-defined]
    Art,
    BaseModelParameters,
    BaseModelState,
    ChildModelParameters,
    ChildModelState,
    Children,
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
    initial_state = _initialise_state(0, state)
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
    state_prev = _initialise_state(i - 1, state)
    state_now = _initialise_state(i, state)

    project_single_year_cpp(i, params, state_prev, state_now)


def run_leapfrog(
    parameters: dict[str, np.ndarray], *, sim_years=None, hts_per_year: int = 10
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
        "hc_nosocomial",
        "hc1_cd4_dist",
        "hc_cd4_transition",
        "hc1_cd4_mort",
        "hc2_cd4_mort",
        "hc1_cd4_prog",
        "hc2_cd4_prog",
        "ctx_effect",
        "ctx_val",
        "hc_art_elig_age",
        "hc_art_elig_cd4",
        "hc_art_mort_rr",
        "hc1_art_mort",
        "hc2_art_mort",
        "hc_art_isperc",
        "hc_art_val",
        "hc_art_init_dist",
        "adult_cd4_dist",
        "fert_mult_by_age",
        "fert_mult_off_art",
        "fert_mult_on_art",
        "total_fertility_rate",
        "local_adj_factor",
        "PMTCT",
        "vertical_transmission_rate",
        "PMTCT_transmission_rate",
        "PMTCT_dropout",
        "PMTCT_input_is_percent",
        "breastfeeding_duration_art",
        "breastfeeding_duration_no_art",
        "mat_hiv_births",
        "mat_prev_input",
        "prop_lt200",
        "prop_gte350",
        "incrate",
        "ctx_val_is_percent",
        "hc_art_is_age_spec",
        "hc_age_coarse",
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

    children = Children(
        params["hc_nosocomial"],
        params["hc1_cd4_dist"],
        params["hc_cd4_transition"],
        params["hc1_cd4_mort"],
        params["hc2_cd4_mort"],
        params["hc1_cd4_prog"],
        params["hc2_cd4_prog"],
        params["ctx_effect"].item(),
        params["ctx_val"][0, :],
        params["hc_art_elig_age"],
        params["hc_art_elig_cd4"],
        params["hc_art_mort_rr"],
        params["hc1_art_mort"],
        params["hc2_art_mort"],
        params["hc_art_isperc"],
        params["hc_art_val"],
        params["hc_art_init_dist"],
        params["adult_cd4_dist"],
        params["fert_mult_by_age"],
        params["fert_mult_off_art"],
        params["fert_mult_on_art"],
        params["total_fertility_rate"],
        params["local_adj_factor"].item(),
        params["PMTCT"],
        params["vertical_transmission_rate"],
        params["PMTCT_transmission_rate"],
        params["PMTCT_dropout"],
        params["PMTCT_input_is_percent"],
        params["breastfeeding_duration_art"],
        params["breastfeeding_duration_no_art"],
        params["mat_hiv_births"],
        params["mat_prev_input"],
        params["prop_lt200"],
        params["prop_gte350"],
        params["incrate"],
        params["ctx_val_is_percent"][0, :],
        params["hc_art_is_age_spec"],
        params["hc_age_coarse"],
    )
    child_model_params = ChildModelParameters(children)
    return Parameters(base_model_params, child_model_params)


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
        p_hiv_deaths=year_state["p_hiv_deaths"],
    )
    child_model_state = ChildModelState(
        hc1_hiv_pop=year_state["hc1_hiv_pop"],
        hc2_hiv_pop=year_state["hc2_hiv_pop"],
        hc1_art_pop=year_state["hc1_art_pop"],
        hc2_art_pop=year_state["hc2_art_pop"],
        hc1_noart_aids_deaths=year_state["hc1_noart_aids_deaths"],
        hc2_noart_aids_deaths=year_state["hc2_noart_aids_deaths"],
        hc1_art_aids_deaths=year_state["hc1_art_aids_deaths"],
        hc2_art_aids_deaths=year_state["hc2_art_aids_deaths"],
        hiv_births=year_state["hiv_births"],
        hc_art_total=year_state["hc_art_total"],
        hc_art_init=year_state["hc_art_init"],
        hc_art_need_init=year_state["hc_art_need_init"],
        ctx_need=year_state["ctx_need"],
        ctx_mean=year_state["ctx_mean"],
    )
    return State(base_model_state, child_model_state)
