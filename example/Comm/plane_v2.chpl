/* This is an example code that attempts to
   quantify how fast one can do the gather/scatters
   necessary for a 3D distributed FFT
*/
use BlockDist;
use Time;

// The grid size
config const Ng=256;
assert((Ng%numLocales)==0,"numLocales should divide Ng");
const SlabSize = Ng/numLocales;

// Define the domain for the grid. Note that the Ng+2 on the last
// dimension is how FFTW handles the real to complex transforms and how
// it packs in the complex data.
const Space = {0.. #Ng, 0.. #Ng, 0.. #Ng};
// Set up the slab decomposition of the array
const targetLocales = reshape(Locales, {0.. #numLocales, 0..0, 0..0});
const D : domain(3) dmapped Block(boundingBox=Space,
                                  targetLocales=targetLocales)=Space;


var arr : [D] complex;
arr = 1.0;

var tt : Timer;
tt.clear();

tt.start();
// Loop over locales
coforall loc in Locales {
  on loc {
    var locD = {here.id*SlabSize.. #SlabSize};
    // Here is the pencil loop
    //forall (j,k) in locD with (var myPencil : [0.. #Ng] complex) {
    var myPlane : [0.. #Ng, 0..0, 0.. #Ng] complex;
    for j in locD {
      // Pull down the data
      myPlane = arr[{0.. #Ng,j..j,0.. #Ng}];

      // Do the 1D FFTs here
      myPlane *= 2;

      // Push back the pencil here
      arr[{0.. #Ng,j..j,0.. #Ng}] = myPlane;
    }
  }
}
tt.stop();


// We should do some validation at the end here.
const diff = max reduce abs(arr)-2.0;

writef("Numlocales=%3i Ng=%6i  ",numLocales, Ng);
writef("  Max diff : %6.2er ",diff);
writef("  Elapsed time  : %8.3r\n", tt.elapsed());

