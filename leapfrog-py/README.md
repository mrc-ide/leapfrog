# leapfrog-py

[![PyPI - Version](https://img.shields.io/pypi/v/leapfrog-py.svg)](https://pypi.org/project/leapfrog-py)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/leapfrog-py.svg)](https://pypi.org/project/leapfrog-py)

-----

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Development](#development)
- [License](#license)

## Installation

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

I use [hatch](https://hatch.pypa.io/1.9/) to manage the project, but this should work without it if you prefer.

If you're using an IDE you might need to set the Python interpreter to the one in the hatch virtual environment. See instructions for [VSCode here](https://hatch.pypa.io/latest/how-to/integrate/vscode/).

We use `Eigen` for linear algebra library, install this on linux

```console
sudo apt-get install libeigen3-dev
```

on windows download from https://eigen.tuxfamily.org/index.php?title=Main_Page, extract the archive and build and install it using cmake. From the extracted dir:

```console
mkdir build
cd build
cmake ..
cmake --build . --target install
```

### Building, installing and running tests

Use hatch

```console
hatch shell
hatch run install
hatch run test
hatch run cov
hatch run install_and_test
hatch run lint:fmt
hatch run lint:style
hatch run lint:typing
hatch run lint:all
```

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
