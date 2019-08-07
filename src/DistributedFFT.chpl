/* Documentation for DistributedFFT */
prototype module DistributedFFT {

  use BlockDist;
  use Barriers;
  use ReplicatedVar;
  use RangeChunk;
  use FFTW;
  use FFTW.C_FFTW;
  require "npFFTW.h";

  extern proc isNullPlan(plan : fftw_plan) : c_int;

  init_FFTW_MT();

  proc deinit() {
    cleanup_threads();
    cleanup();
  }

  config const numFFTWThreads = here.numPUs();
  config const debugTranspose = false;

  var fftw_planner_lock$ : [rcDomain] sync bool;

  enum FFTtype {DFT, R2R};

  /* fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany, */
  /*                              fftw_complex *in, const int *inembed, */
  /*                              int istride, int idist, */
  /*                              fftw_complex *out, const int *onembed, */
  /*                              int ostride, int odist, */
  /*                              int sign, unsigned flags); */
  record FFTWplan {
    var plan : fftw_plan;

    // Mimic the advanced interface 
    proc init(param ftType : FFTtype, numThreads : integral, args ...?k) {
      fftw_planner_lock$(1).writeEF(true);
      fftw_plan_with_nthreads(numThreads:c_int);
      select ftType {
          when FFTtype.DFT do plan = fftw_plan_many_dft((...args));
          when FFTtype.R2R do plan = fftw_plan_many_r2r((...args));
        }
      fftw_planner_lock$(1).readFE();
    }

    proc deinit() {
      destroy_plan(plan);
    }

    proc execute() {
      FFTW.execute(plan);
    }

    proc isValid : bool {
      return isNullPlan(plan)==0;
    }
  }

  /* Convenience constructor for grids */
  proc newSlabDom(dom: domain) where isRectangularDom(dom) {
    if dom.rank < 2 then compilerError("The domain must be at least 2D");
    if ((dom.dim(1).size)%numLocales !=0) then halt("numLocales must divide the first dimension");
    const targetLocales = reshape(Locales, {0.. #numLocales, 0..0, 0..0});
    return dom dmapped Block(boundingBox=dom, targetLocales=targetLocales);
  }

  proc newSlabDom(sz) where isHomogeneousTupleType(sz.type) {
    var tup : (sz.size)*range;
    for param ii in 1..sz.size do tup(ii) = 0.. #sz(ii);
    return newSlabDom({(...tup)});
  }


  /* Warm up the FFTW planner.
     
     arr is assumed to be a block distributed array.
   */
  proc warmUpPlanner(arr : [?Dom]) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");
    if !isComplexType(arr.eltType) then halt("Code is designed for complex arrays only");

    doFFT_YZ(FFTtype.DFT, arr, true, FFTW_FORWARD, FFTW_MEASURE);
    doFFT_YZ(FFTtype.DFT, arr, true, FFTW_BACKWARD, FFTW_MEASURE);
    if (debugTranspose) {
      doFFT_X_Transposed(FFTtype.DFT, new Transposer(), arr, true, FFTW_FORWARD, FFTW_MEASURE);
      doFFT_X_Transposed(FFTtype.DFT, new Transposer(), arr, true, FFTW_BACKWARD, FFTW_MEASURE);
    } else {
      doFFT_X(FFTtype.DFT, arr, true, FFTW_FORWARD, FFTW_MEASURE);
      doFFT_X(FFTtype.DFT, arr, true, FFTW_BACKWARD, FFTW_MEASURE);
    }
  }


  /* FFT.
   */
  proc doFFT(arr : [?Dom], sign : c_int) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");
    if !isComplexType(arr.eltType) then halt("Code is designed for complex arrays only");

    doFFT_YZ(FFTtype.DFT, arr, false, sign, FFTW_WISDOM_ONLY);
    if (debugTranspose) {
      doFFT_X_Transposed(FFTtype.DFT, new Transposer(), arr, false, sign, FFTW_MEASURE);
    } else {
      doFFT_X(FFTtype.DFT, arr, false, sign, FFTW_MEASURE);
    }


    // End of doFFT
  }


  /* Helper routines. This assumes that the data are block-distributed.

     R2C and C2R transforms are more complicated. Currently not supported.

     args -- are the sign and planner flags.

     Note that if you actually want to do a transform, you should have warmed up
     the planner, and then pass FFTW_WISDOM_ONLY in (since the array is otherwise
     destroyed -- unless you use FFTW_ESTIMATE).
   */
  proc doFFT_YZ(param ftType : FFTtype, arr : [?Dom] ?T, warmUpOnly : bool, args ...?k) {

    // Run on all locales
    coforall loc in Locales {
      on loc {
        // Set up for the yz transforms on each domain.
        const myDom = arr.localSubdomain();

        // Get the x-range to loop over
        const xRange = myDom.dim(1);
        const yRange = myDom.dim(2);
        const zRange = myDom.dim(3);

        /* We have a few options here on how to set this up.

           Since plan creation must be single threaded, and is locked,
           we don't want to generate a plan for every yz plane out here.

           We could use the array-execute feature of FFTW, but that
           might run into issues of alignment.

           So, we try a conservative approach first, by just slicing and
           copying the data into a temporary array, doing the FFT and copying
           back. This is identical to the approach we use in doFFT_X.
        */

        // Write down all the parameters explicitly
        var howmany = 1 : c_int;
        var nn : c_array(c_int, 2);
        nn[0] = yRange.size : c_int;
        nn[1] = zRange.size : c_int;
        var nnp = c_ptrTo(nn[0]);
        var rank = 2 : c_int;
        var stride = 1 : c_int;
        var idist = 1 : c_int;

        if (warmUpOnly) {
          var myplane : [{0..0, yRange, zRange}] T;
          var plan_yz = new FFTWplan(ftType, 1, rank, nnp, howmany, c_ptrTo(myplane),
                                    nnp, stride, idist,
                                    c_ptrTo(myplane), nnp, stride, idist,
                                    (...args));
        } else {
          // Pull down each plane, process and send back
          forall i in xRange with
            // Task private variables
            (var myplane : [{0..0, yRange, zRange}] T,
             var plan_yz = new FFTWplan(ftType,1, rank, nnp, howmany, c_ptrTo(myplane),
                                       nnp, stride, idist,
                                       c_ptrTo(myplane), nnp, stride, idist,
                                       (...args))) 
              {
                // Pull down the data
                myplane = arr[{i..i, yRange, zRange}];

                // Do the yz FFTs here
                plan_yz.execute();

                // Push back the pencil here
                arr[{i..i, yRange, zRange}] = myplane;
              }
        }

        // on loc ends
      }
    }
  }
        
  /* X helper.

     See TODO for question on sign
  */
  proc doFFT_X(param ftType: FFTtype, arr : [?Dom], warmUpOnly : bool, args ...?k) {

    coforall loc in Locales {
      on loc {
        // Set up for the yz transforms on each domain.
        const myDom = arr.localSubdomain();
        const localIndex = myDom.first;

        // Split the y-range. Make this dimension agnostic
        const yChunk = chunk(myDom.dim(2), numLocales, here.id);

        // Set FFTW parameters
        const xRange = Dom.dim(1);
        const zRange = myDom.dim(3);
        var nn = xRange.size : c_int;
        var nnp = c_ptrTo(nn);
        var howmany = myDom.dim(3).size : c_int;
        var rank = 1 : c_int;
        var stride = myDom.dim(3).size : c_int;
        var idist = 1 : c_int;

        if (warmUpOnly) {
          var myplane : [{xRange, 0..0, zRange}] complex;
          var plan_x = new FFTWplan(ftType,1, rank, nnp, howmany, c_ptrTo(myplane),
                                    nnp, stride, idist,
                                    c_ptrTo(myplane), nnp, stride, idist,
                                    (...args));
        } else {


          // Pull down each plane, process and send back
          forall j in yChunk with
            // Task private variables
            // TODO : Type here is complex. Does this make sense always???
            (var myplane : [{xRange, 0..0, zRange}] complex,
             var plan_x = new FFTWplan(ftType,1, rank, nnp, howmany, c_ptrTo(myplane),
                                       nnp, stride, idist,
                                       c_ptrTo(myplane), nnp, stride, idist,
                                       (...args))) 
              {
                // Pull down the data
                myplane = arr[{xRange,j..j,zRange}];

                // Do the 1D FFTs here
                plan_x.execute();

                // Push back the pencil here
                arr[{xRange,j..j,zRange}] = myplane;
              }
        }

        // End of on-loc
      }
    }
  }


  /* X helper with a transpose function

     See TODO for question on sign.

     transposeFunc(in, out, sign : bool)
  */
  proc doFFT_X_Transposed(param ftType: FFTtype, transposeFunc, arr : [?Dom] ?T, warmUpOnly : bool, args ...?k) {

    coforall loc in Locales {
      on loc {
        // Set up for the yz transforms on each domain.
        const myDom = arr.localSubdomain();
        const localIndex = myDom.first;

        // Split the y-range. Make this dimension agnostic
        const yChunk = chunk(myDom.dim(2), numLocales, here.id);

        // Set FFTW parameters
        const xRange = Dom.dim(1);
        const zRange = myDom.dim(3);
        var nn = xRange.size : c_int;
        var nnp = c_ptrTo(nn);
        var howmany = myDom.dim(3).size : c_int;
        var rank = 1 : c_int;
        var stride = 1 : c_int;
        var idist = nn;

        if (warmUpOnly) {
          var myplane : [{zRange, 0..0, xRange}] complex;
          var plan_x = new FFTWplan(ftType,1, rank, nnp, howmany, c_ptrTo(myplane),
                                    nnp, stride, idist,
                                    c_ptrTo(myplane), nnp, stride, idist,
                                    (...args));
        } else {


          // Pull down each plane, process and send back
          forall j in yChunk with
            // Task private variables
            // TODO : Type here is complex. Does this make sense always???
            (var myplane : [{zRange, 0..0, xRange}] complex,
             var myplaneT : [{xRange, 0..0, zRange}] T,
             var plan_x = new FFTWplan(ftType,1, rank, nnp, howmany, c_ptrTo(myplane),
                                       nnp, stride, idist,
                                       c_ptrTo(myplane), nnp, stride, idist,
                                       (...args))) 
              {
                // Pull down the data
                myplaneT = arr[{xRange,j..j,zRange}];
                transposeFunc.fwd(myplaneT, myplane);

                // Do the 1D FFTs here
                plan_x.execute();

                // Push back the pencil here
                transposeFunc.rev(myplane, myplaneT);
                arr[{xRange,j..j,zRange}] = myplaneT;
              }
        }

        // End of on-loc
      }
    }
  }

  record Transposer {
    inline proc fwd(inPlane : [?Dom], outPlane) {
      for (i,j,k) in Dom do outPlane[k,j,i] = inPlane[i,j,k];
    }

    inline proc rev(inPlane : [?Dom], outPlane) {
      for (i,j,k) in Dom do outPlane[k,j,i] = inPlane[i,j,k];
    }

  }


  // End of module
}