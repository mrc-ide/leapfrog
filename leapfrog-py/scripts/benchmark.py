#!/usr/bin/env python

import os
import statistics
import timeit

import numpy as np
from leapfrog_py.leapfrog_py import (
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


def input_file_path(file_name):
    return os.path.join("../inst/standalone_model/data/", file_name)


def parameters():
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
