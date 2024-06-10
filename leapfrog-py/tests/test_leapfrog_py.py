import numpy as np
import os
import pytest
from leapfrog import *

@pytest.fixture
def input_file_path():
    def get_path(file_name):
        current_dir = os.path.dirname(__file__)
        return os.path.join(current_dir, '../../inst/standalone_model/data/', file_name)
    return get_path


def read_standalone_data(file_path):
    with open(file_path, 'r') as file:
        # Read the type
        data_type = file.readline().strip()

        # Map Fortran types to numpy types
        type_map = {
            'double': np.float64,
            'int': np.int32
        }
        if data_type not in type_map:
            raise ValueError(f"Unsupported data type: {data_type}")

        np_type = type_map[data_type]

        # Read the dimensions
        dimensions = tuple(map(int, file.readline().strip().split(',')))

        # Read the data
        data = np.fromfile(file, dtype=np_type, sep=',')

        # Reshape data into Fortran-order array
        array = np.reshape(data, dimensions, order='F')

        return array


def get_time_view(array, time):
    return array[..., time]


def create_state_view(output_state, time, child_state):
    state = {k: get_time_view(v, time) for k, v in output_state.items()}
    return State(BaseModelState(**state), child_state)


def test_can_run_model(input_file_path):
    adult_on_art = read_standalone_data(input_file_path("adults_on_art"))
    adults_on_art_is_percent = read_standalone_data(input_file_path("adults_on_art_is_percent"))
    age_specific_fertility_rate = read_standalone_data(input_file_path("age_specific_fertility_rate"))
    art_dropout = read_standalone_data(input_file_path("art_dropout"))
    art_mortality_rate_full = read_standalone_data(input_file_path("art_mortality_rate_full"))
    art_mortality_time_rate_ratio = read_standalone_data(input_file_path("art_mortality_time_rate_ratio"))
    basepop = read_standalone_data(input_file_path("basepop"))
    births_sex_prop = read_standalone_data(input_file_path("births_sex_prop"))
    cd4_initial_distribution_full = read_standalone_data(input_file_path("cd4_initial_distribution_full"))
    cd4_mortality_full = read_standalone_data(input_file_path("cd4_mortality_full"))
    cd4_progression_full = read_standalone_data(input_file_path("cd4_progression_full"))
    idx_hm_elig = read_standalone_data(input_file_path("idx_hm_elig"))
    incidence_age_rate_ratio = read_standalone_data(input_file_path("incidence_age_rate_ratio"))
    incidence_rate = read_standalone_data(input_file_path("incidence_rate"))
    incidence_sex_rate_ratio = read_standalone_data(input_file_path("incidence_sex_rate_ratio"))
    net_migration = read_standalone_data(input_file_path("net_migration"))
    survival_probability = read_standalone_data(input_file_path("survival_probability"))
    h_art_stage_dur = np.array([0.5, 0.5], order='F')

    demog = Demography(basepop, survival_probability, net_migration, age_specific_fertility_rate, births_sex_prop)
    incidence = Incidence(incidence_rate, incidence_age_rate_ratio, incidence_sex_rate_ratio)
    nat_history = NaturalHistory(cd4_mortality_full, cd4_progression_full, cd4_initial_distribution_full, 1)
    art = Art(idx_hm_elig, art_mortality_rate_full, art_mortality_time_rate_ratio, art_dropout, adult_on_art, adults_on_art_is_percent, h_art_stage_dur, 1)
    options = Options(10, 30, 66)
    baseModelParams = BaseModelParameters(options, demog, incidence, nat_history, art)
    childModelParams = BaseModelChildParameters()
    params = Parameters(baseModelParams, childModelParams)

    NS = 2
    pAG = 81
    hAG = 66
    hDS = 7
    hTS = 3
    hAG_span = np.ones(hAG, order='F')
    no_output_years = 2

    p_total_pop = np.zeros((pAG, NS, no_output_years), order='F')
    p_total_pop_natural_deaths = np.zeros((pAG, NS, no_output_years), order='F')
    p_hiv_pop = np.zeros((pAG, NS, no_output_years), order='F')
    p_hiv_pop_natural_deaths = np.zeros((pAG, NS, no_output_years), order='F')
    h_hiv_adult = np.zeros((hDS, hAG, NS, no_output_years), order='F')
    h_art_adult = np.zeros((hTS, hDS, hAG, NS, no_output_years), order='F')
    births = np.zeros((no_output_years), order='F')
    h_hiv_deaths_no_art = np.zeros((hDS, hAG, NS, no_output_years), order='F')
    p_infections = np.zeros((pAG, NS, no_output_years), order='F')
    h_hiv_deaths_art = np.zeros((hTS, hDS, hAG, NS, no_output_years), order='F')
    h_art_initiation = np.zeros((hDS, hAG, NS, no_output_years), order='F')
    p_hiv_deaths = np.zeros((pAG, NS, no_output_years), order='F')
    output_state = {
        "p_total_pop": p_total_pop,
        "p_total_pop_natural_deaths": p_total_pop_natural_deaths,
        "p_hiv_pop": p_hiv_pop,
        "p_hiv_pop_natural_deaths": p_hiv_pop_natural_deaths,
        "h_hiv_adult": h_hiv_adult,
        "h_art_adult": h_art_adult,
        "births": births,
        "h_hiv_deaths_no_art": h_hiv_deaths_no_art,
        "p_infections": p_infections,
        "h_hiv_deaths_art": h_hiv_deaths_art,
        "h_art_initiation": h_art_initiation,
        "p_hiv_deaths": p_hiv_deaths
    }

    child_state = ChildModelState(params)

    state_prev = create_state_view(output_state, 0, child_state)
    print(output_state["p_total_pop"])
    print("Setting initial state")
    set_initial_state(state_prev, params)
    print(output_state["p_total_pop"])
    state_next = create_state_view(output_state, 1, child_state)
    for i in range(1, no_output_years):
        print(f"\n On Step {i}")
        print(state_prev.base.p_total_pop)
        project_single_year(1, params, state_prev, state_next)
        if (i + 1) < no_output_years:
            state_prev = create_state_view(output_state, i, child_state)
            state_next = create_state_view(output_state, i + 1, child_state)

    # print(state_prev.base.p_total_pop)
    print(state_next.base.p_total_pop)

    assert 2 == 2

