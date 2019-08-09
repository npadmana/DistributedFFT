/* Documentation for DistributedFFT */
prototype module DistributedFFT {

  use BlockDist;
  use Barriers;
  use ReplicatedVar;
  use RangeChunk;
  use FFTW;
  use FFTW.C_FFTW;
  use FFT_Timers;
  require "npFFTW.h";


  extern proc isNullPlan(plan : fftw_plan) : c_int;

  proc deinit() {
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
    var tt : TimeTracker;

    // Mimic the advanced interface 
    proc init(param ftType : FFTtype, args ...?k) {
      fftw_planner_lock$(1).writeEF(true);
      //fftw_plan_with_nthreads(numThreads:c_int);
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
      tt.start();
      FFTW.execute(plan);
      tt.stop(TimeStages.Execute);
    }

    proc isValid : bool {
      return isNullPlan(plan)==0;
    }
  }

  /* Convenience constructor for grids */
  proc newSlabDom(dom: domain) where isRectangularDom(dom) {
    // Does this actually work for anything other than 3D? The FFT code does not.
    if dom.rank !=3 then compilerError("The domain must be 3D");
    //if ((dom.dim(1).size)%numLocales !=0) then halt("numLocales must divide the first dimension");
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

    // Note the UNALIGNED flag -- since we do each yz plane separately, make no assumptions
    // about the alignment.
    doFFT_YZ(FFTtype.DFT, arr, true, FFTW_FORWARD, FFTW_MEASURE | FFTW_UNALIGNED);
    doFFT_YZ(FFTtype.DFT, arr, true, FFTW_BACKWARD, FFTW_MEASURE | FFTW_UNALIGNED);
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

    var tt = new TimeTracker();

    tt.start();
    // Must call with WISDOM_ONLY and UNALIGNED
    // WISDOM -- to prevent array from being overwritten; you must call the warmup routine
    // UNALIGNED -- since we do each plane separately, make no alignment assumtions
    doFFT_YZ(FFTtype.DFT, arr, false, sign, FFTW_WISDOM_ONLY | FFTW_UNALIGNED);
    tt.stop(TimeStages.YZ);

    tt.start();
    if (debugTranspose) {
      doFFT_X_Transposed(FFTtype.DFT, new Transposer(), arr, false, sign, FFTW_MEASURE);
    } else {
      doFFT_X(FFTtype.DFT, arr, false, sign, FFTW_MEASURE);
    }
    tt.stop(TimeStages.X);


    // End of doFFT
  }

  /* Helper routines. This assumes that the data are block-distributed.

     R2C and C2R transforms are more complicated. Currently not supported.

     args -- are the sign and planner flags.

     We do each plane separately copying it back and forth.
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

           We use the array-execute interface here and just FFT each plane separately.

           Note that the calling code sets the FFTW_UNALIGNED flag to ensure no alignment mismatches.
        */

        // Write down all the parameters explicitly
        var howmany = 1 : c_int;
        var nn : c_array(c_int, 2);
        nn[0] = yRange.size : c_int;
        nn[1] = zRange.size : c_int;
        var nnp = c_ptrTo(nn[0]);
        var rank = 2 : c_int;
        var stride = 1 : c_int;
        var idist = 0 : c_int;
        var csize : int;
        select T {
            when real do csize=8;
            when complex do csize=16;
            otherwise halt("Unknown type "+T:string);
        }
        csize *= nn[0]*nn[1];
        var arr0 = c_ptrTo(arr.localAccess[myDom.first]);
        var plan_yz = new FFTWplan(ftType, rank, nnp, howmany, arr0,
                                   nnp, stride, idist,
                                   arr0, nnp, stride, idist,
                                   (...args));

        if (!plan_yz.isValid) then
          halt("Error! Plan generation failed! Did you call the warmup routine?");
        if (!warmUpOnly) {
          forall i in xRange {
            var elt = c_ptrTo(arr.localAccess[i,yRange.first, zRange.first]);
            select ftType {
                when FFTtype.DFT do fftw_execute_dft(plan_yz.plan, elt, elt);
                when FFTtype.R2R do fftw_execute_r2r(plan_yz.plan, elt, elt);
              }
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
          var plan_x = new FFTWplan(ftType, rank, nnp, howmany, c_ptrTo(myplane),
                                    nnp, stride, idist,
                                    c_ptrTo(myplane), nnp, stride, idist,
                                    (...args));
        } else {


          // Pull down each plane, process and send back
          forall j in yChunk with
            // Task private variables
            // TODO : Type here is complex. Does this make sense always???
            (var myplane : [{xRange, 0..0, zRange}] complex,
             var plan_x = new FFTWplan(ftType, rank, nnp, howmany, c_ptrTo(myplane),
                                       nnp, stride, idist,
                                       c_ptrTo(myplane), nnp, stride, idist,
                                       (...args)),
             var tt = new TimeTracker()) 
              {
                // Pull down the data
                tt.start();
                myplane = arr[{xRange,j..j,zRange}];
                tt.stop(TimeStages.Comms);

                // Do the 1D FFTs here
                plan_x.execute();

                // Push back the pencil here
                tt.start();
                arr[{xRange,j..j,zRange}] = myplane;
                tt.stop(TimeStages.Comms);
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
          var plan_x = new FFTWplan(ftType, rank, nnp, howmany, c_ptrTo(myplane),
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
             var plan_x = new FFTWplan(ftType, rank, nnp, howmany, c_ptrTo(myplane),
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

  module FFT_Timers {
    use Time;
    // Time the various FFT steps.
    config const timeTrackFFT=false;

    enum TimeStages {X, YZ, Execute, Comms};
    const stageDomain = {TimeStages.X..TimeStages.Comms};
    private var _globalTimeArr : [stageDomain] atomic real;

    resetTimers();

    proc deinit() {
      if timeTrackFFT then printTimers();
    }

    proc resetTimers() {
      for stage in _globalTimeArr do stage.write(0.0);
    }

    proc printTimers() {
      writeln("--------- Timings ---------------");
      for stage in TimeStages {
        writef("Time for %s : %10.2dr\n",stage:string, _globalTimeArr[stage].read());
      }
      writeln("---------------------------------");
    }


    record TimeTracker {
      var tt : Timer();
      var arr : [stageDomain] real;

      proc deinit() {
        if timeTrackFFT {
          for istage in arr.domain {
            _globalTimeArr[istage].add(arr[istage]);
          }
        }
      }

      proc start() {
        if timeTrackFFT {
          tt.clear(); tt.start();
        }
      }

      proc stop(stage) {
        if timeTrackFFT {
          tt.stop();
          arr[stage] += tt.elapsed();
        }
      }
    }

  }


  // End of module
}