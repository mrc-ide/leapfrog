# leapfrog-py

[![PyPI - Version](https://img.shields.io/pypi/v/leapfrog-py.svg)](https://pypi.org/project/leapfrog-py)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/leapfrog-py.svg)](https://pypi.org/project/leapfrog-py)

-----

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Development](#development)
- [License](#license)

## Installation from PyPI

```console
pip install leapfrog-py
```

## Usage

You can use `leapfrog-py` to run a leapfrog simulation in two ways. Run for multiple years using:

```python
from leapfrog_py import run_leapfrog

run_leapfrog(parameters)
```

This will run from 1970 to 2030 inclusive by default with 10 HIV time steps per year. You can also run for a single year using:

```python
from leapfrog_py import set_initial_state, project_single_year

set_initial_state(parameters, state)
for i in range(1, 61):
    project_single_year(i, parameters, state)
```

Parameters and state are both dictionaries of numpy arrays.

## Development

### Prerequisites

This project uses [scikit-build-core](https://github.com/scikit-build/scikit-build-core) to build the C++ project. You'll need a recent version of CMake (>3.15) and Python (>3.7).

I use [uv](https://docs.astral.sh/uv/) to manage the project, but this should work without it if you prefer.

We use `Eigen` for linear algebra library, this will be installed automatically by CMake when this package is built. You can install manually:

On Linux

```console
sudo apt-get install libeigen3-dev
```

On Windows download from https://eigen.tuxfamily.org/index.php?title=Main_Page, extract the archive and build and install it using cmake. From the extracted dir:

```console
mkdir build
cd build
cmake ..
cmake --build . --target install
```

### Building, installing and running tests

* Sync the virtual env `uv sync`
* Run tests `uv run pytest`
* Build wheels `uv build` (currently not working as relying on files from parent)
* Run coverage `uv run pytest --cov --cov-config=pyproject.toml`
* Run lint style `uvx ruff check .` or `uvx black --check --diff .`
* Run lint auto-format `uvx ruff check --fix` or `uvx black .`
* Run lint typing `uv run --group check mypy --install-types --non-interactive src tests`

uv will rebuild automatically if you make a change to the C++ code in the `src` directory, but you will need to force it to reinstall if you updated the C++ library. To force recompilation, pass `--reinstall-package leapfrog-py`. If you want to see verbose output pass `--verbose`. These can be passed and used with any command.

To build with pipx

```console
pipx run build
```

Or simply (in a virtualenv, ideally)

```console
pip install .
```

## License

`leapfrog-py` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
