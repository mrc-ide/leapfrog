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
