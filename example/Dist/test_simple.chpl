use DistributedFFT;
use Random;
use Time;

config const Ng=4;

writef("Running with Ng=%i on %i locales....\n",Ng,numLocales);

const Dom = newSlabDom((Ng,Ng,Ng));
var arr, arr_save : [Dom] complex;
var arr_local : [{0.. #Ng, 0.. #Ng, 0.. #Ng}] complex;

// Local plan
var local_plan = plan_dft(arr_local, FFTW_FORWARD, FFTW_ESTIMATE);

var timeit : Timer;

timeit.clear(); timeit.start();
warmUpPlanner(arr);
timeit.stop();
writef("Time to warm up planner.... %r\n",timeit.elapsed());

fillRandom(arr_local, 11);
fillRandom(arr, 11);

timeit.clear(); timeit.start();
execute(local_plan);
timeit.stop();
writef("Time to locally execute FFT = %r\n",timeit.elapsed());

timeit.clear(); timeit.start();
doFFT(arr, FFTW_FORWARD);
timeit.stop();
writef("Time to run the FFT = %r\n",timeit.elapsed());

arr_save = arr_local;
const maxerr = max reduce abs(arr - arr_save);
writef("The difference between the local and full transform is = %er\n",maxerr);

