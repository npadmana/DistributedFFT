use DistributedFFT;
use Random;
use Time;

config const Ng=64;

const Dom = newSlabDom((Ng,Ng,Ng));
var arr, arrT : [Dom] complex;

var timeit : Timer;

timeit.clear(); timeit.start();
doFFT_Transposed(FFTtype.DFT, arr, arrT, FFTW_FORWARD);
doFFT_Transposed(FFTtype.DFT, arrT, arr, FFTW_BACKWARD);
timeit.stop();
const planTime=timeit.elapsed();

arr = 1.0;

timeit.clear(); timeit.start();
doFFT_Transposed(FFTtype.DFT, arr, arrT, FFTW_FORWARD);
timeit.stop();
const runFtime=timeit.elapsed();

arrT[0,0,0] -= Ng**3;
const maxerr1 = max reduce abs(arrT);
arrT[0,0,0] = Ng**3;

timeit.clear(); timeit.start();
doFFT_Transposed(FFTtype.DFT, arrT, arr, FFTW_BACKWARD);
timeit.stop();
const runBtime=timeit.elapsed();

arr *= 1.0/Ng**3;
const maxerr2 = max reduce abs(arr - 1.0);


writef("numLocales=%2i Ng=%4i diff=%8.2er,%8.2er planTime=%8.2r runTime(F)=%8.2r runTime(B)=%8.2r\n",
       numLocales,Ng, maxerr1, maxerr2, planTime, runFtime, runBtime);