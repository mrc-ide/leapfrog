## Benchmarking

The benchmark script can be used to run the different variants of the leapfrog model via R with timings. Initially this script was used to compare new leapfrog (frogger) and old leapfrog (leaptfrog) timings. Since we have replaced the old version of leapfrog with this repo, we'll no longer be able to run them separately. This document has some timings from two machines for posterity.

On laptop with AMD Ryzen 7 7840U

```
Benchmarking coarse stratified model
# A tibble: 2 × 13
  expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
  <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
1 old          2ms 2.34ms      433.    6.46MB        0    10     0     23.1ms
2 new          2.69ms 2.78ms      360.    1.06MB        0    10     0     27.8ms
# ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
Benchmarking full model
# A tibble: 2 × 13
  expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
  <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
1 old          10.1ms 11.4ms      88.9    4.32MB     9.88     9     1      101ms
2 new          15.8ms 16.6ms      59.5    4.33MB     6.61     9     1      151ms
# ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
Benchmarking child model
# A tibble: 1 × 13
  expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
  <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
1 new          18.9ms 19.1ms      50.6    5.97MB     21.7     7     3      138ms
# ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
```

On my desktop with i7-8700

```
Benchmarking coarse stratified model
# A tibble: 2 × 13
  expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
  <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
1 old          2.37ms  2.7ms      376.    6.35MB        0    10     0     26.6ms
2 new          4.26ms 4.28ms      231.    1.06MB        0    10     0     43.2ms
# ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
Benchmarking full model
# A tibble: 2 × 13
  expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
  <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
1 old          13.3ms 14.9ms      69.9    4.33MB     7.76     9     1      129ms
2 new          21.5ms 21.6ms      45.5    4.33MB     5.06     9     1      198ms
# ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
Benchmarking child model
# A tibble: 1 × 13
  expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
  <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
1 new          24.1ms 24.1ms      41.4    5.58MB     4.60     9     1      217ms
# ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
```
