use DistributedFFT;
use Random;
use Time;
use CommDiagnostics;

config const Ng=4;

writef("Running with Ng=%i on %i locales....\n",Ng,numLocales);

// Use different sizes in x,y,z to catch bugs in the code.
const Dom  = newSlabDom((Ng,   2*Ng, 4*Ng));
const DomT = newSlabDom((2*Ng, Ng,   4*Ng)); // Transposed domain
var arr, arrTT : [Dom] complex; // arrTT will store the inverse of the forward transform
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
doFFT_Transposed(FFTtype.DFT, arrT, arrTT, FFTW_BACKWARD);
timeit.stop();
stopCommDiagnostics();
writef("Time to run the FFT = %r\n",timeit.elapsed());
writef("Note : this includes planning time!\n");

// Validate (both original and final results)
var errBefore, errAfter, errAfter2 : real;
coforall loc in Locales with (max reduce errBefore,
                              max reduce errAfter,
                              max reduce errAfter2) do on loc {
  // Validate original array
  {
    var dom = arr.localSubdomain();
    var tmp = locarr; // This is a little wasteful, but these are small arrays
    const norm = 1.0/(8 * Ng**3);
    var err, err2 : real;
    forall ijk in dom with (max reduce err,
                            max reduce err2) {
      err = abs(tmp[ijk]-arr.localAccess[ijk]);
      err2 = abs(tmp[ijk]-arrTT.localAccess[ijk]*norm);
    }
    errBefore = err;
    errAfter2 = err2;
  }

  // Validate transformed array
  {
    var dom = arrT.localSubdomain();
    var tmp = locarrT;
    var err : real;
    forall (j,i,k) in dom with (max reduce err) {
      err = abs(tmp[i,j,k]-arrT.localAccess[j,i,k]); // Transposed
    }
    errAfter = err;
  }
}

writef("Difference between input array before/after transform = %er\n",errBefore);
writef("Difference between local and distributed transformed arrays = %er\n",errAfter);
writef("Difference between original and transformed-back-and-forth arrays = %er\n",errAfter2);



// Write out comm diagnostics
var comms = getCommDiagnostics();
for ii in 0.. #numLocales {
  writeln("Loc ", ii," ", comms[ii]);
}

