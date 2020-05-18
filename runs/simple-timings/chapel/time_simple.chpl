use DistributedFFT;
use Random;
use Time;

config const Ng=64;
config const loops=1;

const Dom = newSlabDom((Ng,Ng,Ng));
var arr_orig, arr, arrT : [Dom] complex;

var timeit : Timer;

/**** Warmup ****/
timeit.clear(); timeit.start();
doFFT_Transposed(FFTtype.DFT, arr, arrT, FFTW_FORWARD);
doFFT_Transposed(FFTtype.DFT, arrT, arr, FFTW_BACKWARD);
timeit.stop();
const warmupTime=timeit.elapsed();
writeln("Initial warm up complete");

fillRandom(arr_orig, seed=12321);
arr = arr_orig;
writeln("Random array initialized...");

/******** Timing begins here *********************/
var runFtime, runBtime : real;
const norm = 1.0/(Ng:real ** 3);

runFtime=0.0;
runBtime=0.0;

for iloop in 1..loops {
  timeit.clear(); timeit.start();
  doFFT_Transposed(FFTtype.DFT, arr, arrT, FFTW_FORWARD);
  timeit.stop();
  runFtime+=timeit.elapsed();

  timeit.clear(); timeit.start();
  doFFT_Transposed(FFTtype.DFT, arrT, arr, FFTW_BACKWARD);
  timeit.stop();
  runBtime+=timeit.elapsed();

  // Normalize
  forall ijk in Dom {
    arr.localAccess[ijk] *= norm;
  }
}
runFtime /= loops;
runBtime /= loops;

/********** Validate *******************/
var maxdiff=1.0e30;
forall ijk in Dom with (min reduce maxdiff) {
  maxdiff = abs(arr_orig.localAccess[ijk] - arr.localAccess[ijk]);
}


writef("numLocales=%2i Ng=%4i diff=%8.2er planTime=%8.2r runTime(F)=%8.2r runTime(B)=%8.2r\n",
       numLocales,Ng, maxdiff, warmupTime, runFtime, runBtime);