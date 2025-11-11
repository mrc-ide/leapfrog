## Standalone use of leapfrogs model fit

### Important

This standalone version of the model exists only as a proof of concept and
should not actually be run this way. We
anticipate that non-R users will want to build leapfrog as a library (e.g.
using `nanobind`) and call directly from a
program not from this CLI application. This standalone code just demonstrates
one way that we can build leapfrog without
any dependency on R.

## Building with cmake

### Prerequisites

You will need hdf5 installed as this is the serialisation format we use.
To install on linux run

```
sudo apt-get install libhdf5-dev
```

### CMake build

To build with cmake run

```
cmake -B build .
cmake --build build
```

## Running the fit

See usage with

```
./build/simulate_model --help
```

Run as

```
./build/simulate_model 61 ../leapfrog/tests/testthat/testdata/adult_parms_full.h5 output
```

Where

* 1st arg is number of sim years, 61 for 1970:2030 inclusive
* 2nd arg is the path to the input data, relative to this dir
* 3rd arg is the path to the dir where output should be saved

the output will be a hdf5 file named `output.h5` within the specified directory

