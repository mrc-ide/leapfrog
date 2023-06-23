## Standalone use of froggers model fit

### Important

This standalone version of the model exists only as a proof of concept and should not actually be run this way. We
anticipate that non-R users will want to build frogger as a library (e.g. using `pybind11`) and call directly from a
program not from this CLI application. This standalone code just demonstrates one way that we can build frogger without
any dependency on R.

Similarly, our serialisation format here is just the simplest thing that could work and not something we expect we will
support long term; it's just a starting point for a discussion.

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

## Profiling

There are 3 profile runners here, I've only checked these work on linux. Note it might be worth commenting out the step
in the standalone runner which saves the output to disk. This takes a long time and won't be run in the actual model fit
process. Use the valgrind profiler to measure CPU instructions, the google profiler and perf profiler runs a sampling
profiler so should be better for understanding time bottlenecks. I'm not sure really which is the best to use here! Perf
profiler can give a nice view of the callchain which can be expanded and explored.

1. `profile_valgrind` which runs the [callgrind valgrind](https://valgrind.org/docs/manual/cl-manual.html) tool to
   create a profile and then open it with `kcachegrind`.
   You need to have valgrind and kcachegrind installed, then run the script `./profile_valgrind`
2. `profile_gperftoold` uses google [gperftools](https://github.com/gperftools/gperftools) to run the fit
   and [google-pprof](https://github.com/google/pprof) to visualise it.
    * Install `sudo apt-get install libgoogle-perftools-dev`
    * You need `go` to install google-pprof
      ```
      sudo apt-get install golang-go
      go install github.com/google/pprof@latest
      ```
3. `profile_perf`, you will need to install perf
    * `apt-get install linux-tools-common linux-tools-generic linux-tools-`uname -r``
    * Set the perf paranoia level with `sudo sh -c 'echo 1 >/proc/sys/kernel/perf_event_paranoid'`
      see https://superuser.com/questions/980632/run-perf-without-root-rights
    * Run the script as sudo `sudo ./profile_perf`
