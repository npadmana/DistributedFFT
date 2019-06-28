/* Documentation for DistributedFFT */
prototype module DistributedFFT {

  use BlockDist;
  use Barriers;
  use ReplicatedVar;
  use RangeChunk;
  use FFTW;
  use FFTW.C_FFTW;

  init_FFTW_MT();

  proc deinit() {
    cleanup_threads();
    cleanup();
  }

  config const numFFTWThreads = here.numPUs();

  var fftw_planner_lock$ : [rcDomain] sync bool;

  /* fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany, */
  /*                              fftw_complex *in, const int *inembed, */
  /*                              int istride, int idist, */
  /*                              fftw_complex *out, const int *onembed, */
  /*                              int ostride, int odist, */
  /*                              int sign, unsigned flags); */
  record FFTWplan {
    var plan : fftw_plan;

    // Mimic the advanced interface 
    proc init(numThreads : integral, args ...?k) {
      fftw_planner_lock$(1).writeEF(true);
      fftw_plan_with_nthreads(numThreads:c_int);
      plan = fftw_plan_many_dft((...args));
      fftw_planner_lock$(1).readFE();
    }

    proc deinit() {
      destroy_plan(plan);
    }

    proc execute() {
      FFTW.execute(plan);
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
  proc warmUpPlanner(arr : [?Dom] ) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");
    if !isComplexType(arr.eltType) then halt("Code is designed for complex arrays only");

    coforall loc in Locales {
      on loc {
        // Set up for the yz transforms on each domain.
        const myDom = arr.localSubdomain();
        const localIndex = myDom.first;

        // Write down all the parameters explicitly
        var howmany : c_int = myDom.dim(1).size : c_int;
        var nn : c_array(c_int, 2);
        var rank = 2 : c_int;
        var stride = 1 : c_int;
        nn[0] = myDom.dim(2).size : c_int;
        nn[1] = myDom.dim(3).size : c_int;
        var idist = (nn[0]*nn[1]):c_int;
        var nnp = c_ptrTo(nn[0]);
        /* fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany, */
        /*                              fftw_complex *in, const int *inembed, */
        /*                              int istride, int idist, */
        /*                              fftw_complex *out, const int *onembed, */
        /*                              int ostride, int odist, */
        /*                              int sign, unsigned flags); */
        var plan1 = new FFTWplan(numFFTWThreads, rank, nnp, howmany, c_ptrTo(arr[localIndex]),
                                 nnp, stride, idist,
                                 c_ptrTo(arr[localIndex]), nnp, stride, idist,
                                 FFTW_FORWARD, FFTW_MEASURE);
        var plan2 = new FFTWplan(numFFTWThreads, rank, nnp, howmany, c_ptrTo(arr[localIndex]),
                                 nnp, stride, idist,
                                 c_ptrTo(arr[localIndex]), nnp, stride, idist,
                                 FFTW_BACKWARD, FFTW_MEASURE);

        // Warm up the x transpose versions next
        const xRange = Dom.dim(1);
        var plane : [{xRange, myDom.dim(3)}] complex;
        nn[0] = xRange.size : c_int;
        howmany = myDom.dim(3).size : c_int;
        rank = 1 : c_int;
        stride = myDom.dim(3).size : c_int;
        idist = 1 : c_int;
        var plan3 = new FFTWplan(1, rank, nnp, howmany, c_ptrTo(plane),
                                 nnp, stride, idist,
                                 c_ptrTo(plane), nnp, stride, idist,
                                 FFTW_FORWARD, FFTW_MEASURE);
        var plan4 = new FFTWplan(1, rank, nnp, howmany, c_ptrTo(plane),
                                 nnp, stride, idist,
                                 c_ptrTo(plane), nnp, stride, idist,
                                 FFTW_BACKWARD, FFTW_MEASURE);
      }
    }
  }


  /* FFT */
  proc doFFT(arr : [?Dom], sign : c_int) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");
    if !isComplexType(arr.eltType) then halt("Code is designed for complex arrays only");

    var wall = new Barrier(numLocales);

    // Run on all locales
    coforall loc in Locales {
      on loc {
        // Set up for the yz transforms on each domain.
        const myDom = arr.localSubdomain();
        const localIndex = myDom.first;

        // Write down all the parameters explicitly
        var howmany : c_int = myDom.dim(1).size : c_int;
        var nn : c_array(c_int, 2);
        var rank = 2 : c_int;
        var stride = 1 : c_int;
        nn[0] = myDom.dim(2).size : c_int;
        nn[1] = myDom.dim(3).size : c_int;
        var idist = (nn[0]*nn[1]):c_int;
        var nnp = c_ptrTo(nn[0]);
        var plan_yz = new FFTWplan(numFFTWThreads, rank, nnp, howmany, c_ptrTo(arr[localIndex]),
                                 nnp, stride, idist,
                                 c_ptrTo(arr[localIndex]), nnp, stride, idist,
                                 sign, FFTW_WISDOM_ONLY);

        // Execute the plan
        plan_yz.execute();

        wall.barrier();

        // Split the y-range. Make this dimension agnostic
        const yChunk = chunk(myDom.dim(2), numLocales, here.id);

        // Set FFTW parameters
        const xRange = Dom.dim(1);
        const zRange = myDom.dim(3);
        nn[0] = xRange.size : c_int;
        howmany = myDom.dim(3).size : c_int;
        rank = 1 : c_int;
        stride = myDom.dim(3).size : c_int;
        idist = 1 : c_int;

        // Pull down each plane, process and send back
        forall j in yChunk with
          // Task private variables
          (var myplane : [{xRange, 0..0, myDom.dim(3)}] complex,
           var plan_x = new FFTWplan(1, rank, nnp, howmany, c_ptrTo(myplane),
                                     nnp, stride, idist,
                                     c_ptrTo(myplane), nnp, stride, idist,
                                     sign, FFTW_MEASURE))
            {
              // Pull down the data
              myplane = arr[{xRange,j..j,zRange}];

              // Do the 1D FFTs here
              plan_x.execute();

              // Push back the pencil here
              arr[{xRange,j..j,zRange}] = myplane;
            }

        // End of on-loc
      }
    }

    // End of doFFT
  }


}