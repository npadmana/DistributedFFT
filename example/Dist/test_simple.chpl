use DistributedFFT;
use Random;
use Time;
use CommDiagnostics;

config const Ng=4;

writef("Running with Ng=%i on %i locales....\n",Ng,numLocales);

// Use different sizes in x,y,z to catch bugs in the code.
const Dom  = newSlabDom((Ng,   2*Ng, 4*Ng));
const DomT = newSlabDom((2*Ng, Ng,   4*Ng)); // Transposed domain
var arr : [Dom] complex;
var arrT : [DomT] complex;
var locarr, locarrT : [{0.. #Ng, 0.. #2*Ng, 0.. #4*Ng}] complex;

// Local plan
var local_plan = plan_dft(locarr, locarrT, FFTW_FORWARD, FFTW_MEASURE);

var timeit : Timer;

// Fill both the distributed and local arrays with the same random numbers.
fillRandom(locarr, 11);
fillRandom(arr, 11);

timeit.clear(); timeit.start();
execute(local_plan);
timeit.stop();
writef("Time to locally execute FFT = %r\n",timeit.elapsed());

startCommDiagnostics();
timeit.clear(); timeit.start();
doFFT_Transposed(FFTtype.DFT, arr, arrT, FFTW_FORWARD);
timeit.stop();
stopCommDiagnostics();
writef("Time to run the FFT = %r\n",timeit.elapsed());
writef("Note : this includes planning time!\n");

// Validate (both original and final results)
var errBefore, errAfter : real;
coforall loc in Locales with (max reduce errBefore, max reduce errAfter) do on loc {
  // Validate original array
  {
    var tmp = locarr; // Copy in here.
    var dom = arr.localSubdomain();
    var err : real;
    forall ijk in dom with (max reduce err) {
      err = abs(tmp[ijk]-arr.localAccess[ijk]);
    }
    errBefore = err;
  }

  // Validate transformed array
  {
    var tmp = locarrT; // Copy in here.
    var dom = arrT.localSubdomain();
    var err : real;
    forall (j,i,k) in dom with (max reduce err) {
      err = abs(tmp[i,j,k]-arrT.localAccess[j,i,k]); // Transposed
    }
    errAfter = err;
  }
}

writef("Difference between input array before/after transform = %er\n",errBefore);
writef("Difference between local and distributed transformed arrays = %er\n",errAfter);



// Write out comm diagnostics
var comms = getCommDiagnostics();
for ii in 0.. #numLocales {
  writeln("Loc ", ii," ", comms[ii]);
}

