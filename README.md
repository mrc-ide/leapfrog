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

## TODO
* Restructuring the model code to identify more common code
   * There are examples like general demographic projection and hiv population demographic projection which are running similar processes like ageing, non HIV mortality, migration. We should be able to write a function for e.g. ageing which we can run on each of our population matrices. Even for the HIV and ART stratified we can add overloaded function to work with higher dimension data
* Add a test that checks that no `double`s are used in `inst/include` dir. We should be using templated `real_type` for TMB
* Add a broad level overview of the algorithm - is there a diagram available?
* Add a test which proves that C++ code can be compiled standalone, one exists in dust that can be used as inspriration
* Tidy data copying a tensor2 to matrix and a tensor n to array utility functions which should be used in C++ R interface "src/frogger.cpp"
* Ensure input data is not copied, we can read from the data that R owns (as long as don't write to it)
* Convert string flags to the model to enums if we need to switch on them in several places. This should make it easier to reason about in C++ world and isolate the string checking to a single place

## License

MIT © Imperial College of Science, Technology and Medicine
