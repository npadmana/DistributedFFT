use DistributedFFT;

const BaseDom={0..2,0..2,0..2};
const Dom = newSlabDom(BaseDom);
var A,At : [Dom] real;
var r2r = [FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00];

forall (i,j,k) in Dom {
  A[i,j,k] = ((i%2) + 3*(j%2) + (k%2)):real;
}

doFFT_Transposed(FFTtype.R2R,A, At,r2r);

writeln(At);

/*

  This is the non-transposed version, which we store here
  for reference

  The only non-zero elements are
  (0,0,0) -> 160.0
  (0,0,2) -> -32.0

  // These two elements get swapped by a transpose.
  (0,2,0) -> -96.0
  (2,0,0) -> -32.0

var expected= reshape([160.0, 0.0, -32.0,
               0.0, 0.0, 0.0,
               -96.0, 0.0, 0.0,

               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,

               -32.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0],BaseDom);
*/

// This is the transposed version
var expected= reshape([160.0, 0.0, -32.0,
               0.0, 0.0, 0.0,
               -32.0, 0.0, 0.0,

               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,

               -96.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0],BaseDom);
At -= expected;
var diff = max reduce abs(At);
if (diff > 1.0e-14) then writeln("Test FAILED");