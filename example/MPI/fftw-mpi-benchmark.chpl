// Modules we need. Note that we use the existing
// Chapel FFTW module for some common routines. We also
// require MPI here. 
use MPI;
use DistributedFFT; 
use SysCTypes;
use ReplicatedVar;
use Time;
require "fftw3-mpi.h";

// Allow the user to define the number of threads used by FFTW.
config const numFFTWThreads=0;


// Initialize the FFTW library
// We have to do this on every locale, since
// locales correspond to MPI ranks.
coforall loc in Locales {
  on loc {

    fftw_init_threads();
    fftw_mpi_init();

    const nth = if numFFTWThreads > 0 then numFFTWThreads : c_int
      else here.maxTaskPar : c_int;
    fftw_plan_with_nthreads(nth);
  }
}


/* Now that we have everything initialized, we can
   start setting up the necessary arrays.

   FFTW assumes a slab-decomposition of the array, which for 3D grids
   implies that we distribute along the first dimension of the array.
   We use the Block distribution to distribute the array.

   For this example, let us assume that the grid dimension is `{Ng, Ng, Ng}`.

   FFTW suggests using the `fftw_mpi_local_size*` functions to figure
   out how much memory is necessary for FFTW. However, if `Ng` is divisible
   by the number of locales, then we don't need to do anything special, so we
   use that here.
*/

// The grid size
config const Ng=256;
assert((Ng%numLocales)==0,"numLocales should divide Ng");

// Define the domain for the grid. Note that the Ng+2 on the last
// dimension is how FFTW handles the real to complex transforms and how
// it packs in the complex data.
const D = newSlabDom((Ng,Ng,Ng));

// Define the arrays
var A : [D] complex;
A = 1.0;


/*
 Now construct the plan for the forward transform.
*/
var plan : [rcDomain] fftw_plan;

var timer : Timer;
timer.clear(); timer.start();
coforall loc in Locales {
  on loc {
    Barrier(CHPL_COMM_WORLD);

    // Get a pointer to the local part of the array
    const localIndex = (A.localSubdomain()).first;
    var Aptr = c_ptrTo(A.localAccess[localIndex]);

    // This plan is local to the locale
    plan(1) = fftw_mpi_plan_dft_3d(Ng : c_ptrdiff,
                                   Ng : c_ptrdiff,
                                   Ng : c_ptrdiff,
                                   Aptr, Aptr,
                                   CHPL_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
    Barrier(CHPL_COMM_WORLD);
  }
}
timer.stop();
const planTime = timer.elapsed();


timer.clear(); timer.start();
coforall loc in Locales {
  on loc {
    Barrier(CHPL_COMM_WORLD);

    // Run the FFT
    execute(plan(1));
    Barrier(CHPL_COMM_WORLD);
  }
}
timer.stop();
const runTime = timer.elapsed();


timer.clear(); timer.start();
coforall loc in Locales {
  on loc {
    Barrier(CHPL_COMM_WORLD);
    // Clean up
    destroy_plan(plan(1));
    Barrier(CHPL_COMM_WORLD);
  }
}
timer.stop();
const destroyTime = timer.elapsed();

writef("numLocales=%2i Ng=%4i planTime=%8.2r runTime=%8.2r destroyTime=%8.2r\n",
       numLocales,Ng,planTime, runTime, destroyTime);


// Here is where we clean up.
// Again, we must run this on all locales.
coforall loc in Locales {
  on loc {
    fftw_mpi_cleanup();
    fftw_cleanup_threads();
  }
}


// These are the external C declarations
extern const FFTW_MPI_TRANSPOSED_IN : c_uint;
extern const FFTW_MPI_TRANSPOSED_OUT : c_uint;
extern proc fftw_mpi_init();
extern proc fftw_mpi_cleanup();
extern proc fftw_mpi_plan_dft_3d(n0 : c_ptrdiff, n1 : c_ptrdiff, n2 : c_ptrdiff,
                                 inarr : c_ptr(complex) , outarr : c_ptr(complex),
                                 comm : MPI_Comm, sign : c_int, flags : c_uint) : fftw_plan;
