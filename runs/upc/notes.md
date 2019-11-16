Note that this is what I've been able to figure out
based on reading the README and bits of the code.

```
module swap PrgEnv-intel PrgEnv-cray
module load cray-fftw
```

The UPC nomenclature doesn't exactly map on to the NPB.FT class nomenclature

NPB Class    <--> UPC Class <--> Size
D            <--> DD        <--> 2048x1024x1024
E            <--> DD8       <--> 4096x2048x2048
F            <--> DD64      <--> 8192x4096x4096

These are defined in params.h

```
make upc-bench CLASS=DD
make upc-bench CLASS=DD8
make upc-bench CLASS=DD64
```

Quoting from the Running portion of the code below


RUNNING
=======

The benchmark must be run using 8192 UPC threads.

There are two command line parameters (nx, ny) to represent the 
number of upc threads in the X and Y direction of a thread 
partition grid (nx * ny = total number of upc threads). 
Users are allowed to change the values of nx and ny.

For example, on on the Cray at NERSC, the command line will be:
aprun -n 8192 ./ft-2d-upc.fftw3.DD16 64 128

where nx=64 and ny=128 (nx*ny = 8192)

Depending on your system type and configuration, it may be necessary to change
the environment variable XT_SYMMETRIC_HEAP_SIZE 
to define the shared heap size.  A sample run script is provided.
