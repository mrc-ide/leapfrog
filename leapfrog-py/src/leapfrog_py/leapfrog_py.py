import numpy as np

from leapfrog import (  # type: ignore[attr-defined]
    Art,
    ChildModelParameters,
    ChildModelState,
    Children,
    DemographicProjectionParameters,
    DemographicProjectionState,
    Demography,
    HivSimulationParameters,
    HivSimulationState,
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
        "demography_base_pop",
        "demography_survival_probability",
        "demography_net_migration",
        "demography_age_specific_fertility_rate",
        "demography_births_sex_prop",
        "incidence_total_rate",
        "incidence_relative_risk_age",
        "incidence_relative_risk_sex",
        "naturalhistory_cd4_mortality",
        "naturalhistory_cd4_progression",
        "naturalhistory_cd4_initial_distribution",
        "naturalhistory_scale_cd4_mortality",
        "art_idx_hm_elig",
        "art_mortality",
        "art_mortality_time_rate_ratio",
        "art_dropout_recover_cd4",
        "art_dropout_rate",
        "art_adults_on_art",
        "art_adults_on_art_is_percent",
        "art_h_art_stage_dur",
        "art_initiation_mortality_weight",
        "children_hc_nosocomial",
        "children_hc1_cd4_dist",
        "children_hc_cd4_transition",
        "children_hc1_cd4_mort",
        "children_hc2_cd4_mort",
        "children_hc1_cd4_prog",
        "children_hc2_cd4_prog",
        "children_ctx_effect",
        "children_ctx_val",
        "children_hc_art_elig_age",
        "children_hc_art_elig_cd4",
        "children_hc_art_mort_rr",
        "children_hc1_art_mort",
        "children_hc2_art_mort",
        "children_hc_art_isperc",
        "children_hc_art_val",
        "children_hc_art_init_dist",
        "children_adult_cd4_dist",
        "children_fert_mult_by_age",
        "children_fert_mult_off_art",
        "children_fert_mult_on_art",
        "children_total_fertility_rate",
        "children_local_adj_factor",
        "children_PMTCT",
        "children_vertical_transmission_rate",
        "children_PMTCT_transmission_rate",
        "children_PMTCT_dropout",
        "children_PMTCT_input_is_percent",
        "children_breastfeeding_duration_art",
        "children_breastfeeding_duration_no_art",
        "children_mat_hiv_births",
        "children_mat_prev_input",
        "children_prop_lt200",
        "children_prop_gte350",
        "children_incrate",
        "children_ctx_val_is_percent",
        "children_hc_art_is_age_spec",
        "children_hc_age_coarse",
        "children_abortion",
        "children_patients_reallocated",
        "children_hc_art_ltfu",
        "children_hc_age_coarse_cd4",
        "children_adult_female_infections",
        "children_adult_female_hivnpop",
        "children_total_births"
    ]
    assert_keys(params, required_input)
    assert_fortran_order(params)

    demography = Demography(
        params["demography_base_pop"],
        params["demography_survival_probability"],
        params["demography_net_migration"],
        params["demography_age_specific_fertility_rate"],
        params["demography_births_sex_prop"],
    )

    demographic_model_params = DemographicProjectionParameters(
        demography
    )

    incidence = Incidence(
        params["incidence_total_rate"],
        params["incidence_relative_risk_age"],
        params["incidence_relative_risk_sex"],
    )
    nat_history = NaturalHistory(
        params["naturalhistory_cd4_mortality"],
        params["naturalhistory_cd4_progression"],
        params["naturalhistory_cd4_initial_distribution"],
        params["naturalhistory_scale_cd4_mortality"].item(),
    )
    art = Art(
        params["art_idx_hm_elig"],
        params["art_mortality"],
        params["art_mortality_time_rate_ratio"],
        params["art_dropout_recover_cd4"].item(),
        params["art_dropout_rate"],
        params["art_adults_on_art"],
        params["art_adults_on_art_is_percent"],
        params["art_h_art_stage_dur"],
        params["art_initiation_mortality_weight"].item(),
    )

    hiv_simulation_params = HivSimulationParameters(
        incidence, nat_history, art
    )

    children = Children(
        params["children_hc_nosocomial"],
        params["children_hc1_cd4_dist"],
        params["children_hc_cd4_transition"],
        params["children_hc1_cd4_mort"],
        params["children_hc2_cd4_mort"],
        params["children_hc1_cd4_prog"],
        params["children_hc2_cd4_prog"],
        params["children_ctx_effect"].item(),
        params["children_ctx_val"][0, :],
        params["children_hc_art_elig_age"],
        params["children_hc_art_elig_cd4"],
        params["children_hc_art_mort_rr"],
        params["children_hc1_art_mort"],
        params["children_hc2_art_mort"],
        params["children_hc_art_isperc"],
        params["children_hc_art_val"],
        params["children_hc_art_init_dist"],
        params["children_adult_cd4_dist"],
        params["children_fert_mult_by_age"],
        params["children_fert_mult_off_art"],
        params["children_fert_mult_on_art"],
        params["children_total_fertility_rate"],
        params["children_local_adj_factor"].item(),
        params["children_PMTCT"],
        params["children_vertical_transmission_rate"],
        params["children_PMTCT_transmission_rate"],
        params["children_PMTCT_dropout"],
        params["children_PMTCT_input_is_percent"],
        params["children_breastfeeding_duration_art"],
        params["children_breastfeeding_duration_no_art"],
        params["children_mat_hiv_births"],
        params["children_mat_prev_input"],
        params["children_prop_lt200"],
        params["children_prop_gte350"],
        params["children_incrate"],
        params["children_ctx_val_is_percent"][0, :],
        params["children_hc_art_is_age_spec"],
        params["children_hc_age_coarse"],
        params["children_abortion"][0, :, :],
        params["children_patients_reallocated"],
        params["children_hc_art_ltfu"],
        params["children_hc_age_coarse_cd4"],
        params["children_adult_female_infections"],
        params["children_adult_female_hivnpop"],
        params["children_total_births"],
    )
    child_model_params = ChildModelParameters(children)

    # TODO: Move this into C++, 30 is time ART start and 66 is no of HIV age groups which should come from the
    # state space in C++
    options = Options(hts_per_year, 30, 66, 0)
    return Parameters(options, demographic_model_params, hiv_simulation_params, child_model_params)


def _initialise_state(time_step: int, state: dict[str, np.ndarray]) -> State:
    year_state = {k: np.atleast_1d(v[..., time_step]) for k, v in state.items()}
    demographic_projection_stae = DemographicProjectionState(
        p_total_pop=year_state["p_total_pop"],
        births=year_state["births"],
        p_total_pop_natural_deaths=year_state["p_total_pop_natural_deaths"]
    )
    hiv_simulation_state = HivSimulationState(
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
        hc_art_init=year_state["hc_art_init"],
        hc_art_need_init=year_state["hc_art_need_init"],
        ctx_need=year_state["ctx_need"],
        ctx_mean=year_state["ctx_mean"],
    )
    return State(demographic_projection_stae, hiv_simulation_state, child_model_state)
