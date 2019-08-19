Performance comparisons of Chapel vs Reference MPI NPB-FT

Scalability Hardware (Cray-XC):
 - 36-core (72 HT), 128 GB RAM
   - dual 18-core (36 HT) "Broadwell" 2.1 GHz processors

Software:
 - CLE 7.0.UP01
 - intel 19.0.4.243
 - cray-fftw 3.3.8.4.3
 - cray-mpich 7.7.10.2
 - chpl 1.20.0 pre(2487b396fc)

NPB Versions:
 - https://www.nas.nasa.gov/publications/npb.html (NPB3.4-MPI)
 - https://github.com/npadmana/DistributedFFT (ft_transposed.chpl @1ffdfd2) [1]

Note that MPI is only using 32 cores per node)

Dir structure:
 - graphs/      -- contains .gpi files and generated .pdfs
   - data/      -- contains data collated into .dat files
   - raw-data/  -- contains output from runs
     - print.sh -- grep time/perf (modify for size and time/perf)
     - run.sh   -- script I used to run problem sizes


Results in table form:
======================

Size D:
-------

| nodes | Ref MPI | Chapel  |
| ----- | ------- | ------- |
|   1   | 262.25s | 158.37s |
|   2   | 156.53s |  85.18s |
|   4   |  93.80s |  57.99s |
|   8   |  51.06s |  36.10s |
|  16   |  25.49s |  20.59s |
|  32   |  14.11s |  14.07s |
|  64   |  10.45s |   9.39s |
| 128   |   5.25s |   7.23s |
| 256   |   2.82s |   6.09s |


Size E:
-------

| nodes | Ref MPI | Chapel  |
| ----- | ------- | ------- |
|   8   | 510.72s | 291.06s |
|  16   | 240.13s | 160.80s |
|  32   | 115.82s |  96.38s |
|  64   |  59.01s |  63.01s |
| 128   |  51.82s |  42.62s |
| 256   |  25.77s |  30.37s |


Size F:
-------

| nodes | Ref MPI | Chapel  |
| ----- | ------- | ------- |
|  32   | OOM     | OOM     |
|  64   | 974.02s | OOM     |
| 128   | 457.08s | 315.17s |
| 256   | 266.70s | 194.09s |
| 512   | 133.74s | 138.23s |


[1]: https://github.com/npadmana/DistributedFFT/blob/1ffdfd2/example/NPB-FT/ft_transposed.chpl
