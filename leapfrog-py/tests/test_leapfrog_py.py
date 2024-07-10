import os
import subprocess

import numpy as np
import pytest
from leapfrog_py import (
    project_single_year,
    run_leapfrog,
    set_initial_state,
)


@pytest.fixture
def parameters():
    current_dir = os.path.dirname(__file__)
    test_data_dir = os.path.join(
        current_dir, "../../inst/standalone_model/data/child_data"
    )
    test_data_files = [
        f
        for f in os.listdir(test_data_dir)
        if os.path.isfile(os.path.join(test_data_dir, f))
    ]

    parameters = {
        f: read_tensor(os.path.join(test_data_dir, f)) for f in test_data_files
    }
    parameters["art_h_art_stage_dur"] = np.array([0.5, 0.5], order="F")

    # Serialized data is from R, so if we have any data which references an index it will be off by 1
    # so we need to convert it here
    parameters["art_idx_hm_elig"] -= 1
    parameters["children_hc_art_elig_cd4"] -= 1

    return parameters


@pytest.fixture
def state():
    NS = 2  # noqa: N806
    pAG = 81  # noqa: N806
    hAG = 66  # noqa: N806
    hDS = 7  # noqa: N806
    hTS = 3  # noqa: N806
    hc1DS = 7  # noqa: N806
    hc2DS = 6  # noqa: N806
    hc1AG = 5  # noqa: N806
    hc2AG = 10  # noqa: N806
    hcTT = 4  # noqa: N806
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
        "hc1_hiv_pop": np.zeros(
            (hc1DS, hcTT, hc1AG, NS, no_output_years), order="F"
        ),
        "hc2_hiv_pop": np.zeros(
            (hc2DS, hcTT, hc2AG, NS, no_output_years), order="F"
        ),
        "hc1_art_pop": np.zeros(
            (hTS, hc1DS, hc1AG, NS, no_output_years), order="F"
        ),
        "hc2_art_pop": np.zeros(
            (hTS, hc2DS, hc2AG, NS, no_output_years), order="F"
        ),
        "hc1_noart_aids_deaths": np.zeros(
            (hc1DS, hcTT, hc1AG, NS, no_output_years), order="F"
        ),
        "hc2_noart_aids_deaths": np.zeros(
            (hc2DS, hcTT, hc2AG, NS, no_output_years), order="F"
        ),
        "hc1_art_aids_deaths": np.zeros(
            (hTS, hc1DS, hc1AG, NS, no_output_years), order="F"
        ),
        "hc2_art_aids_deaths": np.zeros(
            (hTS, hc2DS, hc2AG, NS, no_output_years), order="F"
        ),
        "hiv_births": np.zeros(no_output_years, order="F"),
        "hc_art_total": np.zeros((4, no_output_years), order="F"),
        "hc_art_init": np.zeros((4, no_output_years), order="F"),
        "hc_art_need_init": np.zeros(
            (hc1DS, hcTT, 15, NS, no_output_years), order="F"
        ),
        "ctx_need": np.zeros(no_output_years, order="F"),
        "ctx_mean": np.zeros(no_output_years, order="F"),
    }


def read_tensor(file_path):
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

    assert np.all(
        state["p_total_pop"][:, :, 0] == parameters["demography_base_pop"]
    )


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


def test_single_year_and_full_fit_agree(parameters, state):
    full_fit = run_leapfrog(parameters)

    set_initial_state(parameters, state)
    for i in range(1, 61):
        project_single_year(i, parameters, state)

    for key in full_fit:
        assert np.all(
            full_fit[key] == state[key]
        ), f"Different result for output {key}"


def test_year_bounds_are_checked(parameters, state):
    with pytest.raises(
        ValueError,
        match="Year must be greater than 0, 0th year is the initial state. See 'set_initial_state'.",
    ):
        project_single_year(0, parameters, state)


def test_can_run_full_model(parameters):
    out = run_leapfrog(parameters)

    # TODO: this is duplicating tests from R, ideally we could compare the output to the R output

    # There is HIV population after age 15
    assert np.all(out["p_hiv_pop"][15:, :, 60] > 0)

    # Natural deaths start at index 17 as no deaths in first HIV population
    # projection as they are calculated from the no of HIV +ve in previous year
    assert np.all(out["p_hiv_pop_natural_deaths"][16:, :, 60] != 0)

    # Some of older ages can be 0 p_infections, so check the middle chunk
    assert np.all(out["p_infections"][15:69, :, 60] > 0)

    assert np.all(out["h_hiv_adult"][:, :, :, 60] != 0)

    # Outputs cannot be negative
    for key in out:
        assert np.all(out[key] >= 0)


def test_can_run_full_model_for_specified_years(parameters):
    out = run_leapfrog(parameters, sim_years=np.arange(1970, 1976))

    assert out["p_hiv_pop"].shape == (81, 2, 6)


def test_produces_same_output_as_r(tmp_path, parameters):
    out = run_leapfrog(parameters)

    current_dir = os.path.dirname(__file__)
    r_scripts_dir = os.path.join(current_dir, "../../scripts")
    subprocess.run(
        [  # noqa: S603 S607
            "Rscript",
            os.path.join(r_scripts_dir, "run_model.R"),
            "--output-dir",
            tmp_path,
        ],
        check=True,
        capture_output=True,
    )

    r_output_files = [
        f
        for f in os.listdir(tmp_path)
        if os.path.isfile(os.path.join(tmp_path, f))
    ]
    assert len(r_output_files) == len(out)

    r_output = {
        f: read_tensor(os.path.join(tmp_path, f)) for f in r_output_files
    }

    assert out.keys() == r_output.keys()

    for key in out:
        # Check they are close up to some tolerance
        # precision will be lost when the R script serializes the data
        assert np.all(
            np.isclose(out[key], r_output[key], atol=1e-5)
        ), f"Different result for output {key}"
