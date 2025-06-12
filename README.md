# frogger

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/mrc-ide/frogger/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/mrc-ide/frogger/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Leapfrog is a multistate population projection model for demographic and
HIV epidemic estimation with interfaces in R, Python, C++ and Delphi.

The name *leapfrog* is in honor of
[Professor](https://iussp.org/en/basia-zaba-1949-2018) Basia
[Zaba](https://translate.google.co.uk/?sl=pl&tl=en&text=Zaba&op=translate).

## Installation

### R

Please install from our
[r-universe](https://mrc-ide.r-universe.dev/builds):

``` r
install.packages(
  "frogger",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
```

You can install the development version of frogger from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("mrc-ide/frogger", subdir = "R")
```

### Python

Note the Python interface is a work in progress

Install from pypi

```
pip install leapfrog-py
```

## Simulation model

The simulation model is implemented in a header-only C++ library located
in [`r-package/int/include/frogger.hpp`](r-package/int/include/frogger.hpp). This location
allows the C++ code to be imported in other R packages via specifying
`LinkingTo: leapfrog` in the `DESCRIPTION` file.

> [!IMPORTANT]
> We use C++20 for this package. Please make sure you have a compiler that is compatible.

See the [R README](r-package/README.md) for details of running the model from R.

## Code design

### Simulation model

The simulation model is implemented as templated C++ code in
`r-package/int/include/frogger.hpp`. This is so the simulation model may be
developed as a standalone C++ library that can be called by other
software without requiring R-specific code features. The code uses
header-only open source libraries to maximize portability.

### Code generation

To change what parameters can be passed in from any interface or the structure of
`Intermediate`, `State` or `OutputState`, please modify json files
[here](./cpp_generation/modelSchemas/).

Then to run code generation follow
[cpp\_generation/README.md](./cpp_generation/README.md)

## Development notes

### Simulation model

  - The model was implemented using *Eigen::Tensor* containers. These
    were preferred for several reasons:
      - Benchmarking found they were slighlty more efficient than
        *boost::multi\_array*.
      - Column-major indexing in the same order as R
      - Other statistical packages (e.g. TMB, Stan) rely heavily on
        *Eigen* so using *Eigen* containers slims the dependencies.

## License

MIT © Imperial College of Science, Technology and Medicine
