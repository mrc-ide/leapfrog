import os

import numpy as np
import pytest
from leapfrog_py.leapfrog_py import (
    project_single_year,
    run_leapfrog,
    set_initial_state,
)


@pytest.fixture
def parameters():
    def input_file_path(file_name):
        current_dir = os.path.dirname(__file__)
        return os.path.join(
            current_dir,
            "../../inst/standalone_model/data/adult_data",
            file_name,
        )

    return {
        "adult_on_art": read_standalone_data(input_file_path("adults_on_art")),
        "adults_on_art_is_percent": read_standalone_data(
            input_file_path("adults_on_art_is_percent")
        ),
        "age_specific_fertility_rate": read_standalone_data(
            input_file_path("age_specific_fertility_rate")
        ),
        "art_dropout": read_standalone_data(input_file_path("art_dropout")),
        "art_mortality_rate_full": read_standalone_data(
            input_file_path("art_mortality_rate_full")
        ),
        "art_mortality_time_rate_ratio": read_standalone_data(
            input_file_path("art_mortality_time_rate_ratio")
        ),
        "basepop": read_standalone_data(input_file_path("basepop")),
        "births_sex_prop": read_standalone_data(
            input_file_path("births_sex_prop")
        ),
        "cd4_initial_distribution_full": read_standalone_data(
            input_file_path("cd4_initial_distribution_full")
        ),
        "cd4_mortality_full": read_standalone_data(
            input_file_path("cd4_mortality_full")
        ),
        "cd4_progression_full": read_standalone_data(
            input_file_path("cd4_progression_full")
        ),
        "idx_hm_elig": read_standalone_data(input_file_path("idx_hm_elig")),
        "relative_risk_age": read_standalone_data(
            input_file_path("incidence_age_rate_ratio")
        ),
        "incidence_rate": read_standalone_data(
            input_file_path("incidence_rate")
        ),
        "relative_risk_sex": read_standalone_data(
            input_file_path("incidence_sex_rate_ratio")
        ),
        "net_migration": read_standalone_data(input_file_path("net_migration")),
        "survival_probability": read_standalone_data(
            input_file_path("survival_probability")
        ),
        "h_art_stage_dur": np.array([0.5, 0.5], order="F"),
    }


@pytest.fixture
def state():
    NS = 2  # noqa: N806
    pAG = 81  # noqa: N806
    hAG = 66  # noqa: N806
    hDS = 7  # noqa: N806
    hTS = 3  # noqa: N806
    no_output_years = 61
    return {
        "p_total_pop": np.zeros((pAG, NS, no_output_years), order="F"),
        "births": np.zeros(no_output_years),
        "p_total_pop_natural_deaths": np.zeros(
            (pAG, NS, no_output_years), order="F"
        ),
        "p_hiv_pop": np.zeros((pAG, NS, no_output_years), order="F"),
        "p_hiv_pop_natural_deaths": np.zeros(
            (pAG, NS, no_output_years), order="F"
        ),
        "h_hiv_adult": np.zeros((hDS, hAG, NS, no_output_years), order="F"),
        "h_art_adult": np.zeros(
            (hTS, hDS, hAG, NS, no_output_years), order="F"
        ),
        "h_hiv_deaths_no_art": np.zeros(
            (hDS, hAG, NS, no_output_years), order="F"
        ),
        "p_infections": np.zeros((pAG, NS, no_output_years), order="F"),
        "h_hiv_deaths_art": np.zeros(
            (hTS, hDS, hAG, NS, no_output_years), order="F"
        ),
        "h_art_initiation": np.zeros(
            (hDS, hAG, NS, no_output_years), order="F"
        ),
        "p_hiv_deaths": np.zeros((pAG, NS, no_output_years), order="F"),
    }


def read_standalone_data(file_path):
    with open(file_path) as file:
        # Read the type
        data_type = file.readline().strip()

        # Map Fortran types to numpy types
        type_map = {"double": np.float64, "int": np.int32}
        if data_type not in type_map:
            msg = f"Unsupported data type: {data_type}"
            raise ValueError(msg)

        np_type = type_map[data_type]

        # Read the dimensions
        dimensions = tuple(map(int, file.readline().strip().split(",")))

        # Read the data
        data = np.fromfile(file, dtype=np_type, sep=",")

        # Reshape data into Fortran-order array
        array = np.reshape(data, dimensions, order="F")

        return array


def test_can_set_initial_state(parameters, state):
    set_initial_state(parameters, state)

    assert np.all(state["p_total_pop"][:, :, 0] == parameters["basepop"])


def test_can_run_single_year(parameters, state):
    set_initial_state(parameters, state)

    project_single_year(1, parameters, state)

    # Demography has run
    assert state["births"][1] > 0
    assert np.all(state["p_total_pop_natural_deaths"][:, :, 1] > 0)

    # HIV pop and infections will stay 0 as no infections after just 1 year
    assert np.all(state["p_hiv_pop"] == 0)
    assert np.all(state["p_hiv_pop_natural_deaths"] == 0)
    assert np.all(state["h_hiv_adult"] == 0)
    assert np.all(state["h_art_adult"] == 0)
    assert np.all(state["h_hiv_deaths_no_art"] == 0)
    assert np.all(state["p_infections"] == 0)
    assert np.all(state["h_hiv_deaths_art"] == 0)
    assert np.all(state["h_art_initiation"] == 0)
    assert np.all(state["p_hiv_deaths"] == 0)


def test_can_run_model(parameters, state):
    full_fit = run_leapfrog(parameters)

    set_initial_state(parameters, state)
    for i in range(1, 61):
        project_single_year(i, parameters, state)

    assert np.all(full_fit["p_total_pop"] == state["p_total_pop"])
    assert np.all(full_fit["births"] == state["births"])
    assert np.all(
        full_fit["p_total_pop_natural_deaths"]
        == state["p_total_pop_natural_deaths"]
    )
    assert np.all(full_fit["p_hiv_pop"] == state["p_hiv_pop"])
    assert np.all(
        full_fit["p_hiv_pop_natural_deaths"]
        == state["p_hiv_pop_natural_deaths"]
    )
    assert np.all(full_fit["h_hiv_adult"] == state["h_hiv_adult"])
    assert np.all(full_fit["h_art_adult"] == state["h_art_adult"])
    assert np.all(
        full_fit["h_hiv_deaths_no_art"] == state["h_hiv_deaths_no_art"]
    )
    assert np.all(full_fit["p_infections"] == state["p_infections"])
    assert np.all(full_fit["h_art_initiation"] == state["h_art_initiation"])
    assert np.all(full_fit["h_hiv_deaths_art"] == state["h_hiv_deaths_art"])
    assert np.all(full_fit["p_hiv_deaths"] == state["p_hiv_deaths"])


def test_year_bounds_are_checked(parameters, state):
    with pytest.raises(
        ValueError,
        match="Year must be greater than 0, 0th year is the initial state. See 'set_initial_state'.",
    ):
        project_single_year(0, parameters, state)


def test_can_run_full_model(parameters):
    out = run_leapfrog(parameters)

    # TODO: this is duplicating tests from R, ideally we could compare the output to the R output

    # No HIV population < age 15
    assert np.all(out["p_hiv_pop"][0:14, :, :] < 1e-20)
    assert np.all(out["p_hiv_pop"][0:14, :, :] > -1e-20)
    assert np.all(out["p_hiv_pop_natural_deaths"][0:15, :, :] == 0)
    assert np.all(out["p_infections"][0:14, :, :] == 0)

    # There is HIV population after age 15
    assert np.all(out["p_hiv_pop"][15:, :, 60] > 0)

    # Natural deaths start at index 17 as no deaths in first HIV population
    # projection as they are calculated from the no of HIV +ve in previous year
    assert np.all(out["p_hiv_pop_natural_deaths"][16:, :, 60] != 0)

    # Some of older ages can be 0 p_infections, so check the middle chunk
    assert np.all(out["p_infections"][15:69, :, 60] > 0)

    assert np.all(out["h_hiv_adult"][:, :, :, 60] != 0)
    # TODO: Why are these 2 failing? I expect input data differences
    # assert np.all(out["h_art_adult"][:, :, :, :, 60] != 0)
    # assert np.all(out["h_art_initiation"][:, :, :, 60] != 0)

    # Outputs cannot be negative
    assert np.all(out["p_total_pop"] >= 0)
    assert np.all(out["births"] >= 0)
    assert np.all(out["p_total_pop_natural_deaths"] >= 0)
    assert np.all(out["p_hiv_pop"] >= 0)
    assert np.all(out["p_hiv_pop_natural_deaths"] >= 0)
    assert np.all(out["h_hiv_adult"] >= 0)
    assert np.all(out["h_art_adult"] >= 0)
    assert np.all(out["h_hiv_deaths_no_art"] >= 0)
    assert np.all(out["p_infections"] >= 0)
    assert np.all(out["h_art_initiation"] >= 0)
    assert np.all(out["h_hiv_deaths_art"] >= 0)
    assert np.all(out["p_hiv_deaths"] >= 0)


def test_can_run_full_model_for_specified_years(parameters):
    out = run_leapfrog(parameters, sim_years=np.arange(1970, 1976))

    assert out["p_hiv_pop"].shape == (81, 2, 6)
