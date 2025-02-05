
<!-- README.md is generated from README.Rmd. Please edit that file -->

# frogger

<!-- badges: start -->

[![Project Status: Concept – Minimal or no implementation has been done
yet, or the repository is only intended to be a limited example, demo,
or
proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build
status](https://github.com/mrc-ide/frogger/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/frogger/actions)
[![codecov.io](https://codecov.io/github/mrc-ide/frogger/coverage.svg?branch=main)](https://codecov.io/github/mrc-ide/frogger?branch=main)
<!-- badges: end -->

Leapfrog is a multistate population projection model for demographic and
HIV epidemic estimation.

The name *leapfrog* is in honor of
[Professor](https://iussp.org/en/basia-zaba-1949-2018) Basia
[Zaba](https://translate.google.co.uk/?sl=pl&tl=en&text=Zaba&op=translate).

## Installation

You can install the development version of frogger from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("mrc-ide/frogger")
```

## Simulation model

The simulation model is implemented in a header-only C++ library located
in [`inst/include/frogger.hpp`](inst/include/frogger.hpp). This location
allows the C++ code to be imported in other R packages via specifying
`LinkingTo: leapfrog` in the `DESCRIPTION` file.

The simulation model is callable in R via a wrapper function
`run_model()` created with [Rcpp](https://www.rcpp.org).

You can control how the simulation model is run with the following
arguments:

- `run_hiv_simulation` which is `TRUE` by default. Set to `FALSE` to
  turn off the HIV simulation and run only the demographic projection.
- `hiv_age_stratification` which must be “coarse” or “full”. Coarse is
  run with 5-year age groups and full with single year ages.
- `run_child_model` which is `FALSE` by default. Set to `TRUE` to run
  the child portion of the model.

## Example

The file `pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ`
contains an example Spectrum file constructed from default country data
for Botswana with Spectrum (April 2022).

Prepare model inputs.

``` r
library(frogger)

pjnz <- system.file("pjnz/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ",
                    package = "frogger", mustWork = TRUE)

demp <- prepare_leapfrog_demp(pjnz)
hivp <- prepare_leapfrog_projp(pjnz)
```

Simulate adult ‘full’ age group (single-year age) and ‘coarse’ age group
(collapsed age groups) models from 1970 to 2030 with 10 HIV time steps
per year.

``` r
lsimF <- run_model(demp, hivp, 1970:2030, 10L,
                   hiv_age_stratification = "full", run_child_model = FALSE)
lsimC <- run_model(demp, hivp, 1970:2030, 10L,
                   hiv_age_stratification = "coarse", run_child_model = FALSE)
```

Compare the HIV prevalence age 15-49 years and AIDS deaths 50+ years.
Deaths 50+ years are to show some noticeable divergence between the
`"full"` and `"coarse"` age group simulations.

``` r
prevF <- colSums(lsimF$p_hiv_pop[16:50,,],,2) / colSums(lsimF$p_total_pop[16:50,,],,2)
prevC <- colSums(lsimC$p_hiv_pop[16:50,,],,2) / colSums(lsimC$p_total_pop[16:50,,],,2)

deathsF <- colSums(lsimF$p_hiv_deaths[51:81,,],,2)
deathsC <- colSums(lsimC$p_hiv_deaths[51:81,,],,2)

plot(1970:2030, prevF, type = "l", main = "Prevalence 15-49")
lines(1970:2030, prevC, col = 2)
```

<img src="man/figures/README-sim_prev-1.png" width="100%" />

``` r

plot(1970:2030, deathsF, type = "l", main = "AIDS Deaths 50+ years")
lines(1970:2030, deathsC, col = 2)
```

<img src="man/figures/README-sim_prev-2.png" width="100%" />

## Benchmarking

Install the package and then run the benchmarking script
`./scripts/benchmark`

## lint

Lint R code with `lintr`

``` r
lintr::lint_package()
```

Lint C++ code with [cpplint](https://github.com/cpplint/cpplint)

``` console
cpplint inst/include/*
```

## Code design

### Simulation model

The simulation model is implemented as templated C++ code in
`inst/include/frogger.hpp`. This is so the simulation model may be
developed as a standalone C++ library that can be called by other
software without requiring R-specific code features. The code uses
header-only open source libraries to maximize portability.

### R functions

The file `src/frogger.cpp` contains R wrapper functions for the model
simulation via [Rcpp](http://dirk.eddelbuettel.com/code/rcpp.html) and
[RcppEigen](http://dirk.eddelbuettel.com/code/rcpp.eigen.html).

## Development notes

### Simulation model

- The model was implemented using *Eigen::Tensor* containers. These were
  preferred for several reasons:
  - Benchmarking found they were slighlty more efficient than
    *boost::multi_array*.
  - Column-major indexing in the same order as R
  - Other statistical packages (e.g. TMB, Stan) rely heavily on *Eigen*
    so using *Eigen* containers slims the dependencies.

### TODO

- Restructuring the model code to identify more common code
  - There are examples like general demographic projection and hiv
    population demographic projection which are running similar
    processes like ageing, non HIV mortality, migration. We should be
    able to write a function for e.g. ageing which we can run on each of
    our population matrices. Even for the HIV and ART stratified we can
    add overloaded function to work with higher dimension data
- Add a test that checks that no `double`s are used in `inst/include`
  dir. We should be using templated `real_type` for TMB
- Add a broad level overview of the algorithm - is there a diagram
  available?
- Convert string flags to the model to enums if we need to switch on
  them in several places. This should make it easier to reason about in
  C++ world and isolate the string checking to a single place
- Previously `hiv_negative_pop` was fixed size by having dimensions
  specified by template, how much does this speed up the code? Is there
  a better way to do this?
- Tidy up confusing looping see
  <https://github.com/mrc-ide/frogger/pull/7#discussion_r1217847753>
- Add R casting helpers which return better errors than Rcpp see
  <https://github.com/mrc-ide/frogger/pull/7#discussion_r1217884684>
- Add a helper to do 0 to base 1 conversion and check upper bounds see
  <https://github.com/mrc-ide/frogger/pull/7#discussion_r1217888684>
- Review what we pass as parameters - can some of these be computed in
  the struct ctor?
  e.g. <https://github.com/mrc-ide/frogger/pull/7#discussion_r1217890466>
- Refactor `OutputState` to take a struct of state-space dimensions
  instead of unpacking the subset of parameters we need. See
  <https://github.com/mrc-ide/frogger/pull/12#discussion_r1245170775>

## License

MIT © Imperial College of Science, Technology and Medicine
