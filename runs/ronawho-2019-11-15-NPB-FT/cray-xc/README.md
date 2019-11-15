Performance comparisons of Chapel vs Reference MPI NPB-FT

Scalability Hardware (Cray-XC):
 - 36-core (72 HT), 128 GB RAM
   - dual 18-core (36 HT) "Broadwell" 2.1 GHz processors

Software:
 - CLE 7.0.UP01
 - MPI runs:
   - icc 19.0.4.243
   - cray-fftw 3.3.8.4
   - cray-mpich 7.7.10.2
 - UPC Runs:
   - cce 9.0.2 (classic)
   - cray-fftw 3.3.8.4
 - Chapel Runs:
   - chpl 1.20.0
   - icc 19.0.4.243
   - cray-fftw 3.3.8.4

NPB Versions:
 - https://www.nas.nasa.gov/publications/npb.html (NPB3.4-MPI)
 - http://www.nersc.gov/assets/Trinity--NERSC-8-RFP/Benchmarks/Jan9/UPC-FT.tar
 - https://github.com/npadmana/DistributedFFT (ft_transposed.chpl @9d02ce5) [1]

Note that MPI and UPC are only using 32 cores per node)

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

| nodes | Ref MPI | Ref UPC | Chapel  |
| ----- | ------- | ------- | ------- |
|   1   | 261.57s | -       | 147.20s |
|   2   | 156.37s | -       |  89.49s |
|   4   |  92.73s |  71.26s |  56.54s |
|   8   |  49.99s |  40.79s |  32.63s |
|  16   |  25.32s |  21.54s |  16.99s |
|  32   |  12.77s |  11.93s |  10.40s |
|  64   |  10.48s |   6.31s |   5.74s |
| 128   |   5.19s |   3.66s |   3.12s |
| 256   |   2.59s |   2.94s |   1.88s |
| 512   |   1.26s |   3.30s |   1.02s |


Size E:
-------

| nodes | Ref MPI | Ref UPC | Chapel  |
| ----- | ------- | ------- | ------- |
|   8   | 504.33s | -       | 314.25s |
|  16   | 244.86s | -       | 162.79s |
|  32   | 116.66s | 104.71s |  88.04s |
|  64   |  59.19s |  57.13s |  54.27s |
| 128   |  51.69s |  28.32s |  28.97s |
| 256   |  25.46s |  16.15s |  16.45s |
| 512   |  12.32s |   9.25s |   8.43s |


Size F:
-------
| nodes | Ref MPI | Ref UPC | Chapel  |
| ----- | ------- | ------- | ------- |
|  64   | 977.14s | _       | 403.31s |
| 128   | 452.68s | _       | 242.79s |
| 256   | 261.48s | 152.86s | 137.60s |
| 512   | 133.16s |  72.08s |  70.05s |

[1]: https://github.com/npadmana/DistributedFFT/blob/9d02ce5/example/NPB-FT/ft_transposed.chpl
