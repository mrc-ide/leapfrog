# frogger

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build status](https://github.com/mrc-ide/frogger/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/frogger/actions)
[![codecov.io](https://codecov.io/github/mrc-ide/frogger/coverage.svg?branch=main)](https://codecov.io/github/mrc-ide/frogger?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/frogger/badge)](https://www.codefactor.io/repository/github/mrc-ide/frogger)
<!-- badges: end -->

## Installation

To install `frogger`:

```r
remotes::install_github("mrc-ide/frogger", upgrade = FALSE)
```

## Running

Frogger has different parts of the model that you can turn on or off. At the moment these are
* `hiv_age_stratification` which must be "coarse" or "full". Coarse is run with 5-year age groups and full with single year ages.
* `run_child_model` which is `FALSE` by default but can be set to `TRUE` to run the paediatric portion of the model.

### R

If you are running the model from R just call the R function with the desired arguments e.g. to run with full age stratification and child model turned on.

```R
run_model <- function(data, parameters, sim_years,
                      hts_per_year, output_steps,
                      hiv_age_stratification = "full",
                      run_child_model = TRUE)
```

### C++

If you are running the model from C++, you need to create the `ModelVariant` struct you want to fit and call `run_model` with this template param e.g. to run the child model

```C++
constexpr auto ss = leapfrog::StateSpace<leapfrog::ChildModel>();
ret = leapfrog::run_model<leapfrog::ChildModel>(proj_years, save_steps, params);
```

Other model variants are `leapfrog::BaseModelFullAgeStratification` and `leapfrog::BaseModelCoarseAgeStratification` see `inst/include/model_variants.hpp`

## lint

Lint R code with `lintr`

```R
lintr::lint_package()
```

Lint C++ code with [cpplint](https://github.com/cpplint/cpplint)

```console
cpplint inst/include/*
```

## TODO

* Restructuring the model code to identify more common code
    * There are examples like general demographic projection and hiv population demographic projection which are running
      similar processes like ageing, non HIV mortality, migration. We should be able to write a function for e.g. ageing
      which we can run on each of our population matrices. Even for the HIV and ART stratified we can add overloaded
      function to work with higher dimension data
* Add a test that checks that no `double`s are used in `inst/include` dir. We should be using templated `real_type` for
  TMB
* Add a broad level overview of the algorithm - is there a diagram available?
* Convert string flags to the model to enums if we need to switch on them in several places. This should make it easier to reason about in C++ world and isolate the string checking to a single place
* Update `gender` terminology to sex
* Previously `hiv_negative_pop` was fixed size by having dimensions specified by template, how much does this speed up the code? Is there a better way to do this?
* Tidy up confusing looping see https://github.com/mrc-ide/frogger/pull/7#discussion_r1217847753
* Add R casting helpers which return better errors than Rcpp
  see https://github.com/mrc-ide/frogger/pull/7#discussion_r1217884684
* Add a helper to do 0 to base 1 conversion and check upper bounds
  see https://github.com/mrc-ide/frogger/pull/7#discussion_r1217888684
* Review what we pass as parameters - can some of these be computed in the struct ctor?
  e.g. https://github.com/mrc-ide/frogger/pull/7#discussion_r1217890466
* Refactor `OutputState` to take a struct of state-space dimensions instead of unpacking the subset of parameters we
  need. See https://github.com/mrc-ide/frogger/pull/12#discussion_r1245170775

## License

MIT © Imperial College of Science, Technology and Medicine
