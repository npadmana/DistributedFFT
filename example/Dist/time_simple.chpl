use DistributedFFT;
use Random;
use Time;

config const Ng=64;

const Dom = newSlabDom((Ng,Ng,Ng));
var arr : [Dom] complex;

var timeit : Timer;

timeit.clear(); timeit.start();
warmUpPlanner(arr);
timeit.stop();
const planTime=timeit.elapsed();

arr = 1.0;

timeit.clear(); timeit.start();
doFFT(arr, FFTW_FORWARD);
timeit.stop();
const runFtime=timeit.elapsed();

timeit.clear(); timeit.start();
doFFT(arr, FFTW_BACKWARD);
timeit.stop();
const runBtime=timeit.elapsed();

arr *= 1.0/Ng**3;
const maxerr = max reduce abs(arr - 1.0);


writef("numLocales=%2i Ng=%4i diff=%8.2er planTime=%8.2r runTime(F)=%8.2r runTime(B)=%8.2r\n",
       numLocales,Ng, maxerr, planTime, runFtime, runBtime);