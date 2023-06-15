## Standalone use of froggers model fit

### Using installed frogger & RcppEigen

Install frogger, from an R terminal in this directory run

```r
devtools::install_local("../../")
```

Then run configuration from a shell using

```
./configure
```

which will write out a `Makefile` with the path to your installed copy of frogger's library, then

```
make
```

The resulting program will run the leapfrog model end-to-end saving outputs from the last time point.

### Without installing frogger or RcppEigen

Run configuration script providing the path to `frogger` and `eigen` e.g.

```
./configure ~/projects/frogger /path/to/eigen
```

## Running the fit

See usage with

```
./fit_model --help
```

Run as

```
./fit_model 60 10 data output
```

Where

* 1st arg is number of sim years
* 2nd arg is number of HIV steps within the year
* 3rd arg is the path to the input data, `data` relative to this dir
* 3th arg is the path where output should be saved
