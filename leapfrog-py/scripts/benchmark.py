#!/usr/bin/env python

import os
import statistics
import timeit

import numpy as np

from leapfrog_py import (
    project_single_year,
    run_leapfrog,  # noqa F401
    set_initial_state,
)


def pretty_timeit(stmt, globals, setup="pass", repeat=5, number=1000):
    times = timeit.repeat(
        stmt=stmt, setup=setup, globals=globals, repeat=repeat, number=number
    )

    times_ms = [time * 1000 for time in times]

    min_time = min(times_ms)
    median_time = statistics.median(times_ms)
    max_time = max(times_ms)
    total_time = sum(times_ms)

    output = (
        f"Minimum Time: {min_time:.1f} ms\n"
        f"Median Time: {median_time:.1f} ms\n"
        f"Maximum Time: {max_time:.1f} ms\n"
        f"Number of Iterations: {repeat}\n"
        f"Total Time: {total_time:.1f} ms\n"
    )

    return output


def parameters():
    test_data_dir = "../inst/standalone_model/data/child_data"
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


def fit_by_single(parameters, state):
    set_initial_state(parameters, state)
    for i in range(1, 61):
        project_single_year(i, parameters, state)


params = parameters()
print("Full model benchmark")
print(
    pretty_timeit("run_leapfrog(params)", globals=locals(), number=1, repeat=50)
)
print("Year by year benchmark")
print(
    pretty_timeit(
        "fit_by_single(params, s)",
        globals=locals(),
        setup="s = state()",
        number=1,
        repeat=50,
    )
)
