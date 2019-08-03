/*
  A Chapel implementation of the NPB FT Benchmark, using the
  DistributedFFT routines in this module.

  DistributedFFT implements a slab-decomposed FFT, using FFTW for all the
  local steps, but it does the transpose communication using Chapel.

  The goal of this benchmark is to a) verify that the FFT is indeed correctly
  implemented, and b) to compare its performance to FT benchmark.


*/
use DistributedFFT;
use Time;

// Define the classes
enum NPB {S,W,A,B,C,D,E,F};
config const NPBClass : NPB = NPB.S;

// Constants
const alpha = 1.0e-6;
config const Threshold = 1.0e-12;


// Define the domain
// Note that the definitions in the NPB documentation are
// Fortran-ordered, so we transpose the x and z variables
const (Nz, Ny, Nx) = ProblemSizes(NPBClass);
const Dom = newSlabDom((Nx,Ny,Nz));


// Define arrays
var V, W  : [Dom] complex;
var Twiddle : [Dom] real;

var timeit : Timer;

// Warm up the FFTW planners.
// We don't time this
timeit.clear(); timeit.start();
warmUpPlanner(V);
warmUpPlanner(W);
timeit.stop();
writef("Time to setup FFT plans   : %10.4dr \n",timeit.elapsed());

// Touch the arrays once
timeit.clear(); timeit.start();
initialize_U();
initialize_twiddle();
timeit.stop();
writef("Time to touch arrays once : %10.4dr \n",timeit.elapsed());

// Now do the actual evolution
const refChecksums = ReferenceChecksums(NPBClass);
// Start timers
timeit.clear(); timeit.start();
initialize_U();
initialize_twiddle();

for ii in refChecksums.domain {
  evolve();
  const c1 = checksum();
  const verified = verify(c1, refChecksums[ii]);
  if !verified then halt("Checksum verification failed, halting.....");
  writef("Checksum(%0i) =  %14.12ez\n",ii,c1);
}
timeit.stop();

writeln("Successful completion of NPB FT Benchmark Class ",NPBClass);
const total_time = timeit.elapsed();
const ntotal_f = (Nx*Ny*Nz):real;
const niter = refChecksums.size;
const mflops = 1.0e-6*ntotal_f *
  (14.8157+7.19641*log(ntotal_f)
   +  (5.23518+7.21113*log(ntotal_f))*niter)
  /total_time;
writef("Elapsed time : %10.4dr\n",timeit.elapsed());
writef("MFLOPS : %10.4dr\n",mflops);

////////////////////////////////////////////////////
// Benchmark ends here
///////////////////////////////////////////////////


proc evolve() {
  const FTnorm = 1.0/(Nx*Ny*Nz):real;
  forall ijk in Dom {
    V[ijk] *= Twiddle[ijk];
    W[ijk] = V[ijk];
  }
  doFFT(W, FFTW_BACKWARD);
  W *= FTnorm;
}

proc checksum() : complex {
  var ck : complex = 0.0 + 0.0i;
  forall j in 0..1023 with (+ reduce ck) {
    const q = (5*j)%Nx;
    const r = (3*j)%Ny;
    const s = j%Nz;
    ck += W[q,r,s];
  }
  return ck;
}


proc initialize_twiddle() {
  const N = (Nx, Ny, Nz);
  const halfN = (Nx/2, Ny/2, Nz/2);
  const fac = 4.0*pi**2 * alpha;
  forall ijk in Dom {
    var e = 0.0;
    for param ii in 1..3 {
      const i1 = ijk(ii);
      const x1 = if (i1 >= halfN(ii)) then i1-N(ii) else i1;
      e += x1**2;
    }
    Twiddle[ijk] = exp(-fac*e);
  }
}

// This initializes the U array. Since that is never used,
// we store it in V.
proc initialize_U() {
  use Random;
  fillRandom(V, 314159265, RNG.NPB);
  doFFT(V, FFTW_FORWARD);
}

proc verify(x : complex, y : complex) : bool {
  const diffre = abs((x.re - y.re)/y.re) < Threshold;
  const diffim = abs((x.im - y.im)/y.im) < Threshold;
  return diffre && diffim;
}

///-----------------------------------------------------///
// Problem definitions and checksums go here.

proc ProblemSizes(cc : NPB) : (3*int) {
  select cc {
      when NPB.S do return (64,64,64);
      when NPB.W do return (128, 128, 32);
      when NPB.A do return (256, 256, 128);
      when NPB.B do return (512, 256, 256);
      when NPB.C do return (512, 512, 512);
      when NPB.D do return (2048, 1024, 1024);
      when NPB.E do return (4096, 2048, 2048);
      when NPB.F do return (8192, 4096, 4096);
      otherwise do halt("Unknown class");
  }
}

proc ReferenceChecksums(cc : NPB) {
  select cc {
      when NPB.S do return
                      [(5.546087004964e+02, 4.845363331978e+02): complex,
                       (5.546385409189e+02, 4.865304269511e+02): complex,
                       (5.546148406171e+02, 4.883910722336e+02): complex,
                       (5.545423607415e+02, 4.901273169046e+02): complex,
                       (5.544255039624e+02, 4.917475857993e+02): complex,
                       (5.542683411902e+02, 4.932597244941e+02): complex];
      when NPB.W do return 
                      [(5.673612178944e+02, 5.293246849175e+02): complex,
                       (5.631436885271e+02, 5.282149986629e+02): complex,
                       (5.594024089970e+02, 5.270996558037e+02): complex,
                       (5.560698047020e+02, 5.260027904925e+02): complex,
                       (5.530898991250e+02, 5.249400845633e+02): complex,
                       (5.504159734538e+02, 5.239212247086e+02): complex];
      when NPB.A do return 
                      [(5.046735008193e+02, 5.114047905510e+02): complex,
                       (5.059412319734e+02, 5.098809666433e+02): complex,
                       (5.069376896287e+02, 5.098144042213e+02): complex,
                       (5.077892868474e+02, 5.101336130759e+02): complex,
                       (5.085233095391e+02, 5.104914655194e+02): complex,
                       (5.091487099959e+02, 5.107917842803e+02): complex];
      otherwise do halt("Unknown class");
  }
}



