tenger run: 

 `$ salloc -n 1 -A edu19.summer -t 10:00`

 `$ make test-full`

The results:
```
m = 2050, n = 2050, tEnd = 0.1, num_threads =  1, elapsed = 53.5126
m = 2050, n = 2050, tEnd = 0.1, num_threads =  2, elapsed = 27.0679
m = 2050, n = 2050, tEnd = 0.1, num_threads =  3, elapsed = 23.359
m = 2050, n = 2050, tEnd = 0.1, num_threads =  4, elapsed = 14.5042
m = 2050, n = 2050, tEnd = 0.1, num_threads =  5, elapsed = 14.1487
m = 2050, n = 2050, tEnd = 0.1, num_threads =  6, elapsed = 9.84213
m = 2050, n = 2050, tEnd = 0.1, num_threads =  7, elapsed = 11.2676
m = 2050, n = 2050, tEnd = 0.1, num_threads =  8, elapsed = 7.63202
m = 2050, n = 2050, tEnd = 0.1, num_threads =  9, elapsed = 7.5244
m = 2050, n = 2050, tEnd = 0.1, num_threads = 10, elapsed = 6.66463
m = 2050, n = 2050, tEnd = 0.1, num_threads = 11, elapsed = 6.09392
m = 2050, n = 2050, tEnd = 0.1, num_threads = 12, elapsed = 5.2418
m = 2050, n = 2050, tEnd = 0.1, num_threads = 13, elapsed = 5.47249
m = 2050, n = 2050, tEnd = 0.1, num_threads = 14, elapsed = 4.5869
m = 2050, n = 2050, tEnd = 0.1, num_threads = 15, elapsed = 4.26091
m = 2050, n = 2050, tEnd = 0.1, num_threads = 16, elapsed = 4.23007
m = 2050, n = 2050, tEnd = 0.1, num_threads = 17, elapsed = 3.97215
m = 2050, n = 2050, tEnd = 0.1, num_threads = 18, elapsed = 3.63754
m = 2050, n = 2050, tEnd = 0.1, num_threads = 19, elapsed = 4.02056
m = 2050, n = 2050, tEnd = 0.1, num_threads = 20, elapsed = 3.48301
m = 2050, n = 2050, tEnd = 0.1, num_threads = 21, elapsed = 3.46411
m = 2050, n = 2050, tEnd = 0.1, num_threads = 22, elapsed = 3.12799
m = 2050, n = 2050, tEnd = 0.1, num_threads = 23, elapsed = 4.75923
m = 2050, n = 2050, tEnd = 0.1, num_threads = 24, elapsed = 4.65006
m = 2050, n = 2050, tEnd = 0.1, num_threads = 25, elapsed = 4.57425
m = 2050, n = 2050, tEnd = 0.1, num_threads = 26, elapsed = 4.54472
m = 2050, n = 2050, tEnd = 0.1, num_threads = 27, elapsed = 4.48565
m = 2050, n = 2050, tEnd = 0.1, num_threads = 28, elapsed = 4.40989
m = 2050, n = 2050, tEnd = 0.1, num_threads = 29, elapsed = 4.24804
m = 2050, n = 2050, tEnd = 0.1, num_threads = 30, elapsed = 4.12939
m = 2050, n = 2050, tEnd = 0.1, num_threads = 31, elapsed = 3.97509
m = 2050, n = 2050, tEnd = 0.1, num_threads = 32, elapsed = 3.90429
```

```
$ srun -n 1 orig/shwater2d_naive
Solver took 4.6281 seconds

$ srun -n 1 ./shwater2d 48 2048  
m = 2050, n = 2050, tEnd = 0.1, num_threads = 48, elapsed = 2.79613
```

