# leapfrog-py

## Installation

```bash
pip install leapfrog-py
```

## Usage


### Get parameters

In the following example you will need parameters - inputs to the model. We have provided utilities to easily get test parameters (we use hdf5 file format to store these). If you already have different parameters in the correct shape you can skip this section.

```python
from leapfrog_py import read_h5_file

parameters_adult = read_h5_file("../r-package/tests/testthat/testdata/adult_parms_full.h5")
parameters_child = read_h5_file("../r-package/tests/testthat/testdata/child_parms_full.h5")
```


### Run model

Here are a couple of ways to run the model

```python
from leapfrog_py import run_model

# Default configuration is HivFullAgeStratification and output_years is range(1970, 2031)
run_model(parameters_adult)

# Can give another configuration but make sure the parameters match it
run_model(parameters_child, "ChildModel")

# Can specify which years to output from the model
run_model(parameters_child, "ChildModel", [1970, 1971, 2000])
```


### Run model from initial state

You can also run the model from an initial state. For example, to run the model in two parts

```python
from leapfrog_py import run_model, run_model_from_state, get_time_slice

# Run model until year 2000
ret = run_model(parameters_adult, "HivFullAgeStratification", range(1970, 2001))

run_model_from_state(
    # parameters
    parameters_adult,

    # configuration
    "HivFullAgeStratification",
    
    # get the last time slice of the state returned by run model so we
    # can continue from where the model left off
    get_time_slice(ret, 30),

    # simulation start year, this should be the year the previous
    # argument represents
    2000,

    # years you want the model to output
    range(2001, 2031)
)
```


### Run model for a single year

The previous two functions will return a dictionary with values that are have one more
dimension (time dimension) than the state specified in the [model config](../cpp_generation/modelSchemas/configs/DpConfig.json). This means if you want to carry on from that initial single year state you would have to use the `get_time_slice` utility at index 1. Example:

```python
from leapfrog_py import run_model, run_model_from_state, get_time_slice

ret_single_year = run_model(parameters_adult, "HivFullAgeStratification", [1970])

for year in range(1970, 2031):
    ret_single_year = run_model_from_state(
        parameters_adult,
        "HivFullAgeStratification",
        get_time_slice(ret_first_year, 0),
        year,
        [year + 1]
    )
```

This is tedious for two reasons:
1. Repeatedly getting the time slice for one time dimension (initial one is unavoidable since we use run model and don't have any initial state)
1. Specifying output year as `year + 1`

So to make projecting for a single year more ergonomic we decided to create a separate function
that returns the time slice for a single year:

```python
from leapfrog_py import run_model, run_model_single_year, get_time_slice

ret_single_year = get_time_slice(run_model(parameters_adult, "HivFullAgeStratification", [1970]), 0)

for year in range(1970, 2031):
    ret_single_year = run_model_single_year(
        parameters_adult,
        "HivFullAgeStratification",
        ret_single_year,
        year
    )
```

## Development

This project uses [scikit-build-core](https://github.com/scikit-build/scikit-build-core) with [nanobind](https://github.com/wjakob/nanobind) to build the C++ project. You'll need a recent version of CMake (>3.15) and Python (>3.7).

We use [uv](https://docs.astral.sh/uv/) to manage the project:
* Sync the virtual env `uv sync`
* Run tests `uv run pytest`
* Lint `uvx ruff check .`

uv will rebuild automatically if you make a change to the C++ code however if you need to force recompilation use `--reinstall-package leapfrog-py`.

For local installation run `pip install .`

## License

`leapfrog-py` is distributed under the terms of [MIT](https://spdx.org/licenses/MIT.html) license.
