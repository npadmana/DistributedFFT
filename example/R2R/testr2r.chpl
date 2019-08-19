use DistributedFFT;

const BaseDom={0..2,0..2,0..2};
const Dom = newSlabDom(BaseDom);
var A : [Dom] real;
var r2r = [FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00];

warmUpPlanner(A, r2r);

forall (i,j,k) in Dom {
  A[i,j,k] = ((i%2) + 3*(j%2) + (k%2)):real;
}

doR2R(A, r2r);

writeln(A);

var expected= reshape([160.0, 0.0, -32.0,
               0.0, 0.0, 0.0,
               -96.0, 0.0, 0.0,

               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,

               -32.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0],BaseDom);
A -= expected;
var diff = max reduce abs(A);
if (diff > 1.0e-14) then writeln("Test FAILED");