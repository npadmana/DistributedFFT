[![CI](https://github.com/npadmana/DistributedFFT/actions/workflows/CI.yml/badge.svg)](https://github.com/npadmana/DistributedFFT/actions/workflows/CI.yml)

Documentation for this package may be found [here](https://npadmana.github.io/DistributedFFT/docs/modules/src/DistributedFFT.html).

Branches:
---------

- v0.3: Version of master, as of October 2019. This is before we merged in the elegance 
        branch. For a discussion of the changes + perf tests, see #27.


Tags:
-----
- v0.4: Tag when we merged in the elegance branch (see #27)
- v0.2: Version used in our PAW-ATM 19 submission. Uses single-threaded FFTW w/ Chapel handling 
        multi-threading and some more optimizations. (See #14 for some useful timings)
- v0.1: First pass at working version, with some optimizations. Used multi-threaded FFTW.
