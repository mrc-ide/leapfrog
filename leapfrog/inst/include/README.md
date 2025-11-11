This directory contains header files which can be imported by other R packages using the `LinkingTo:` directive in the `DESCRIPTION` file. See https://r-pkgs.org/src.html?q=linkingto#cpp-import

This directory contains two subdirectories:
* `generated` - this contains automatically generated code and should not be modified by hand, see "Development" vignette for details
* `models` - this directory contains the model implementation code

_Note:_ if this contains a file `leapfrog.h` (the same as the package name), then when building the package `Rcpp::compileAttributes()` will add the line `#include ../include/leapfrog.h` in [`src/RcppExports.cpp`](../../src/RcppExports.cpp). This creates a clash with multiple loading of some Eigen libraries on Windows (only; fine on Linux and Mac environment).

I am unsure why at the moment, but for now avoid `leapfrog.h` until resolved.
