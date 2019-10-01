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

  config const numberOfPlanes=1;

  extern proc isNullPlan(plan : fftw_plan) : c_int;

  proc deinit() {
    cleanup();
  }

  pragma "locale private"
  var fftw_planner_lock$ : chpl_LocalSpinlock;

  enum FFTtype {DFT, R2R};

  /* fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany, */
  /*                              fftw_complex *in, const int *inembed, */
  /*                              int istride, int idist, */
  /*                              fftw_complex *out, const int *onembed, */
  /*                              int ostride, int odist, */
  /*                              int sign, unsigned flags); */
  /* fftw_plan fftw_plan_many_r2r(int rank, const int *n, int howmany, */
  /*                              double *in, const int *inembed, */
  /*                              int istride, int idist, */
  /*                              double *out, const int *onembed, */
  /*                              int ostride, int odist, */
  /*                              const fftw_r2r_kind *kind, unsigned flags); */
  record FFTWplan {
    param ftType : FFTtype;
    var plan : fftw_plan;
    var tt : TimeTracker;

    // Mimic the advanced interface 
    proc init(param ftType1 : FFTtype, args ...?k) {
      ftType = ftType1;
      this.complete();
      fftw_planner_lock$.lock();
      select ftType {
          when FFTtype.DFT do plan = fftw_plan_many_dft((...args));
          when FFTtype.R2R do plan = fftw_plan_many_r2r((...args));
        }
      fftw_planner_lock$.unlock();
    }

    proc deinit() {
      fftw_planner_lock$.lock();
      destroy_plan(plan);
      fftw_planner_lock$.unlock();
    }

    proc execute() {
      tt.start();
      FFTW.execute(plan);
      tt.stop(TimeStages.Execute);
    }

    proc execute(arr1, arr2) {
      select ftType {
          when FFTtype.DFT do fftw_execute_dft(plan, arr1, arr2);
          when FFTtype.R2R do fftw_execute_r2r(plan, arr1, arr2);
        }
    }


    proc isValid : bool {
      return isNullPlan(plan)==0;
    }
  }

  /* Convenience constructor for grids */
  proc newSlabDom(dom: domain) where isRectangularDom(dom) {
    if dom.rank !=3 then compilerError("The domain must be 3D");
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
  proc warmUpPlanner(arr : [?Dom] complex) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");

    // Note the UNALIGNED flag -- since we do each yz plane separately, make no assumptions
    // about the alignment.
    doFFT_YZ(FFTtype.DFT, arr, true, FFTW_FORWARD, FFTW_MEASURE | FFTW_UNALIGNED);
    doFFT_YZ(FFTtype.DFT, arr, true, FFTW_BACKWARD, FFTW_MEASURE | FFTW_UNALIGNED);

    doFFT_X(FFTtype.DFT, arr, true, FFTW_FORWARD, FFTW_MEASURE);
    doFFT_X(FFTtype.DFT, arr, true, FFTW_BACKWARD, FFTW_MEASURE);
  }

  /* Warm up the FFTW planner.
     
     arr is assumed to be a block distributed array.

     Specialize here for R2R transforms.
   */
  proc warmUpPlanner(arr : [?Dom] real, r2rType : [] fftw_r2r_kind) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");
    if r2rType.size != 3 then halt("Need 3 R2R kinds");
    if r2rType.rank != 1 then halt("Expected a 1D array");

    ref r2r = r2rType.reindex(0..2);

    // Note the UNALIGNED flag -- since we do each yz plane separately, make no assumptions
    // about the alignment.
    doFFT_YZ(FFTtype.R2R, arr, true, r2r[1..2], FFTW_MEASURE | FFTW_UNALIGNED);
    doFFT_X(FFTtype.R2R, arr, true, r2r[0..0], FFTW_MEASURE);
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
    doFFT_X(FFTtype.DFT, arr, false, sign, FFTW_MEASURE);
    tt.stop(TimeStages.X);


    // End of doFFT
  }

  /* FFT.

     Stores the FFT in dest transposed (xyz -> yxz).

     Note that both src and dest are overwritten.
   */
  proc doFFT_Transposed(src: [?SrcDom] complex, dest : [?DestDom] complex, sign : c_int) {
    if SrcDom.rank != 3 then halt("Code is designed for 3D arrays only");
    if DestDom.rank != 3 then halt("Code is designed for 3D arrays only");
    if SrcDom.dim(1) != DestDom.dim(2) then halt("Mismatched x-y ranges");
    if SrcDom.dim(2) != DestDom.dim(1) then halt("Mismatched y-x ranges");
    if SrcDom.dim(3) != DestDom.dim(3) then halt("Mismatched z ranges");

    var tt = new TimeTracker();

    tt.start();
    // Must call with WISDOM_ONLY and UNALIGNED
    // WISDOM -- to prevent array from being overwritten; you must call the warmup routine
    // UNALIGNED -- since we do each plane separately, make no alignment assumtions
    doFFT_YZ_Transposed(FFTtype.DFT, src, dest, false, sign, FFTW_MEASURE);
    tt.stop(TimeStages.YZ);
    writef("yz=%dr  ",tt.tt.elapsed());

    tt.start();
    doFFT_X_Transposed(FFTtype.DFT, dest, false, sign, FFTW_MEASURE);
    tt.stop(TimeStages.X);
    writef("x=%dr  ",tt.tt.elapsed());


    // End of doFFT
  }

  /* FFT.

     Stores the FFT in dest transposed (xyz -> yxz).
   */
  proc doFFT_Transposed_Elegant(param ftType : FFTtype,
                                src: [?SrcDom] ?T,
                                dest : [?DestDom] T,
                                signOrKind) {
    // Sanity checks
    if SrcDom.rank != 3 then halt("Code is designed for 3D arrays only");
    if DestDom.rank != 3 then halt("Code is designed for 3D arrays only");
    if SrcDom.dim(1) != DestDom.dim(2) then halt("Mismatched x-y ranges");
    if SrcDom.dim(2) != DestDom.dim(1) then halt("Mismatched y-x ranges");
    if SrcDom.dim(3) != DestDom.dim(3) then halt("Mismatched z ranges");


    var YZBarrier = new Barrier(numLocales);

    coforall loc in Locales {
      on loc {
        const mySrcDom = src.localSubdomain();
        const xSrc = mySrcDom.dim(1);
        const ySrc = mySrcDom.dim(2);
        const zSrc = mySrcDom.dim(3);
        const myLineSize = zSrc.size*c_sizeof(T):int;
        const myDestDom = dest.localSubdomain();
        const yDest = myDestDom.dim(1);
        const xDest = myDestDom.dim(2);

        // Set up FFTW plans
        var plan_x = setup1DPlan(T, ftType, xDest.size, zSrc.size, signOrKind, FFTW_MEASURE);
        var plan_y = setup1DPlan(T, ftType, ySrc.size, zSrc.size, signOrKind, FFTW_MEASURE);
        var plan_z = setup1DPlan(T, ftType, zSrc.size, 1, signOrKind, FFTW_MEASURE);

        // Use this as a temporary work array to
        // avoid rewriting the src array
        var myplane : [{0..0, ySrc, zSrc}] T;

        for ix in xSrc {
          // Copy data over. This is a completely
          // local operation, so use localAccess.
          // Or a series of memcpys for performance.
          [(tmp, iy,iz) in myplane.domain] myplane[0, iy,iz] = src.localAccess[ix,iy,iz];

          // y-transform
          const y0 = ySrc.first;
          forall iz in zSrc with (ref plan_y) {
            var elt = c_ptrTo(myplane[0,y0,iz]);
            plan_y.execute(elt, elt);
          }

          // z-transform
          const z0 = zSrc.first;
          // Offset to reduce collisions
          const offset = (ySrc.size/numLocales)*here.id;
          forall iy in 0.. #ySrc.size with (ref plan_z) {
            const iy1 = (iy + offset)%ySrc.size + y0;
            var elt = c_ptrTo(myplane[0,iy1,z0]);
            plan_z.execute(elt, elt);
            // This is the transpose step
            dest[{iy1..iy1,ix..ix,zSrc}] = myplane[{0..0, iy1..iy1,zSrc}];
          }
        }

        // Wait until all communication is complete
        YZBarrier.barrier();

        // x-transform
        const x0 = xDest.first;
        forall (iy,iz) in {yDest, zSrc} with (ref plan_x) {
          var elt = c_ptrTo(dest.localAccess[iy, x0, iz]);
          plan_x.execute(elt, elt);
        }

        // End of on-loc
      }
    }


    // End of doFFT
  }

  /* R2R
   */
  proc doR2R(arr : [?Dom] real, r2rType : [] fftw_r2r_kind) {
    if arr.rank != 3 then halt("Code is designed for 3D arrays only");
    if r2rType.size != 3 then halt("Need 3 R2R kinds");
    if r2rType.rank != 1 then halt("Expected a 1D array");

    ref r2r = r2rType.reindex(0..2);

    var tt = new TimeTracker();

    tt.start();
    // Must call with WISDOM_ONLY and UNALIGNED
    // WISDOM -- to prevent array from being overwritten; you must call the warmup routine
    // UNALIGNED -- since we do each plane separately, make no alignment assumtions
    // Note the UNALIGNED flag -- since we do each yz plane separately, make no assumptions
    // about the alignment.
    doFFT_YZ(FFTtype.R2R, arr, false, r2r[1..2], FFTW_WISDOM_ONLY | FFTW_UNALIGNED);
    tt.stop(TimeStages.YZ);

    tt.start();
    doFFT_X(FFTtype.R2R, arr, false, r2r[0..0], FFTW_MEASURE);
    tt.stop(TimeStages.X);
    // End of doR2R
  }


  // Set up the YZ transform plan
  proc setupYZPlan(param ftType : FFTtype, yzplane : [?Dom] ?T, signOrKind, flags : c_uint) {
    // Pull signOrKind locally since this may be an array
    // we need to take a pointer to.
    var mySignOrKind = signOrKind;
    var arg0 : _signOrKindType(ftType);
    select ftType {
        when FFTtype.R2R do arg0 = c_ptrTo(mySignOrKind);
        when FFTtype.DFT do arg0 = mySignOrKind;
      }

    // Set up for the yz transforms on each domain.
    const myDom = yzplane.localSubdomain();

    // Get the x-range to loop over
    const yRange = myDom.dim(2);
    const zRange = myDom.dim(3);


    // Write down all the parameters explicitly
    var howmany = 1 : c_int;
    var nn : c_array(c_int, 2);
    nn[0] = yRange.size : c_int;
    nn[1] = zRange.size : c_int;
    var nnp = c_ptrTo(nn[0]);
    var rank = 2 : c_int;
    var stride = 1 : c_int;
    var idist = 0 : c_int;
    var arr0 = c_ptrTo(yzplane.localAccess[myDom.first]);
    return new FFTWplan(ftType, rank, nnp, howmany, arr0,
                        nnp, stride, idist,
                        arr0, nnp, stride, idist,
                        arg0, flags);
  }

  // Set up 1D in-place plans
  proc setup1DPlan(type arrType, param ftType : FFTtype, nx : int, strideIn : int, signOrKind, in flags : c_uint) {
    // Pull signOrKind locally since this may be an array
    // we need to take a pointer to.
    var mySignOrKind = signOrKind;
    var arg0 : _signOrKindType(ftType);
    select ftType {
        when FFTtype.R2R do arg0 = c_ptrTo(mySignOrKind);
        when FFTtype.DFT do arg0 = mySignOrKind;
      }

    // Define a dummy array
    var arr : [0.. #(nx*strideIn)] arrType;

    // Write down all the parameters explicitly
    var howmany = 1 : c_int;
    var nn : c_array(c_int, 1);
    nn[0] = nx : c_int;
    var nnp = c_ptrTo(nn[0]);
    var rank = 1 : c_int;
    var stride = strideIn  : c_int;
    var idist = 0 : c_int;
    var arr0 = c_ptrTo(arr);
    flags = flags | FFTW_UNALIGNED;
    return new FFTWplan(ftType, rank, nnp, howmany, arr0,
                        nnp, stride, idist,
                        arr0, nnp, stride, idist,
                        arg0, flags);
  }

  // Set up many 1D in place plans
  proc setupPlanColumns(type arrType, param ftType : FFTtype, dom : domain(2), numTransforms : int, signOrKind, in flags : c_uint) {
    // Pull signOrKind locally since this may be an array
    // we need to take a pointer to.
    var mySignOrKind = signOrKind;
    var arg0 : _signOrKindType(ftType);
    select ftType {
        when FFTtype.R2R do arg0 = c_ptrTo(mySignOrKind);
        when FFTtype.DFT do arg0 = mySignOrKind;
      }

    // Define a dummy array
    var arr : [dom] arrType;

    // Write down all the parameters explicitly
    var howmany = numTransforms : c_int;
    var nn : c_array(c_int, 1);
    nn[0] = dom.dim(1).size : c_int;
    var nnp = c_ptrTo(nn[0]);
    var rank = 1 : c_int;
    var stride = dom.dim(2).size  : c_int;
    var idist = 1 : c_int;
    var arr0 = c_ptrTo(arr);
    flags = flags | FFTW_UNALIGNED;
    return new FFTWplan(ftType, rank, nnp, howmany, arr0,
                        nnp, stride, idist,
                        arr0, nnp, stride, idist,
                        arg0, flags);
  }


  // Set up 1D out-of-place plans
  proc setup1DPlan(type arrType, param ftType : FFTtype, nx : int, strideIn : int, strideOut : int, signOrKind, in flags : c_uint) {
    // Pull signOrKind locally since this may be an array
    // we need to take a pointer to.
    var mySignOrKind = signOrKind;
    var arg0 : _signOrKindType(ftType);
    select ftType {
        when FFTtype.R2R do arg0 = c_ptrTo(mySignOrKind);
        when FFTtype.DFT do arg0 = mySignOrKind;
      }

    // Define a dummy array
    var arrIn : [0.. #(nx*strideIn)] arrType;
    var arrOut : [0.. #(nx*strideOut)] arrType;

    // Write down all the parameters explicitly
    var howmany = 1 : c_int;
    var nn : c_array(c_int, 1);
    nn[0] = nx : c_int;
    var nnp = c_ptrTo(nn[0]);
    var rank = 1 : c_int;
    var idist = 0 : c_int;
    flags = flags | FFTW_UNALIGNED;
    return new FFTWplan(ftType, rank, nnp, howmany, c_ptrTo(arrIn),
                        nnp, strideIn : c_int, idist,
                        c_ptrTo(arrOut), nnp, strideOut : c_int, idist,
                        arg0, flags);
  }

  /* Helper routines. This assumes that the data are block-distributed.

     R2C and C2R transforms are more complicated. Currently not supported.

     args -- are the sign and planner flags.

     We do each plane separately copying it back and forth.
   */
  proc doFFT_YZ(param ftType : FFTtype, arr : [?Dom] ?T, warmUpOnly : bool, signOrKind, flags : c_uint) {


    // Run on all locales
    coforall loc in Locales {
      on loc {
        // Get the x-range to loop over
        const myDom = arr.localSubdomain();
        const xRange = myDom.dim(1);
        const yRange = myDom.dim(2);
        const zRange = myDom.dim(3);

        var plan_z = setup1DPlan(T, ftType, zRange.size, 1, signOrKind, flags);
        var plan_y = setup1DPlan(T, ftType, yRange.size, zRange.size, signOrKind, flags);

        var tt = new TimeTracker();
        tt.start();
        if (!warmUpOnly) {
          // Do the z transforms
          const z0 = zRange.first;
          forall (ix, iy) in {xRange, yRange} with (ref plan_z) {
            var elt = c_ptrTo(arr.localAccess[ix, iy, z0]);
            plan_z.execute(elt, elt);
          }
          // Do the y transforms
          const y0 = yRange.first;
          forall (ix, iz) in {xRange, zRange} with (ref plan_y) {
            var elt = c_ptrTo(arr.localAccess[ix, y0, iz]);
            plan_y.execute(elt, elt);
          }
        }
        tt.stop(TimeStages.Execute);

        // on loc ends
      }
    }
  }

  proc doFFT_YZ_Transposed(param ftType : FFTtype, arr : [?Dom] ?T, dest : [?DomDest] T, warmUpOnly : bool, signOrKind, flags : c_uint) {


    // Run on all locales
    coforall loc in Locales {
      on loc {
        // Get the x-range to loop over
        const myDom = arr.localSubdomain();
        const xRange = myDom.dim(1);
        const yRange = myDom.dim(2);
        const zRange = myDom.dim(3);
        const myLineSize = zRange.size*c_sizeof(T):int;

        const numTasks = min(here.maxTaskPar, zRange.size);
        const numTransforms = zRange.size/numTasks;
        var plan_y = setupPlanColumns(T, ftType, {yRange, zRange}, numTransforms, signOrKind, flags);
        var plan_y1 = setupPlanColumns(T, ftType, {yRange, zRange}, numTransforms+1, signOrKind, flags);
        var plan_z = setup1DPlan(T, ftType, zRange.size, 1, signOrKind, flags);

        var tt = new TimeTracker();
        tt.start();
        if (!warmUpOnly) {
          // Do the y transforms
          const y0 = yRange.first;

          coforall itask in 0..#numTasks with (ref plan_y, ref plan_y1) {
            const myzRange = chunk(zRange, numTasks, itask);
            for ix in xRange {
              var elt = c_ptrTo(arr.localAccess[ix, y0, myzRange.first]);
              select myzRange.size {
                  when numTransforms do plan_y.execute(elt,elt);
                  when numTransforms+1 do plan_y1.execute(elt,elt);
                  otherwise halt("Bad myzRange size");
                }
            }
          }

          // Do the z transforms
          const z0 = zRange.first;
          const offset = (yRange.size/numLocales)*here.id;
          forall (ix, iy) in {xRange, 0.. #yRange.size} with (ref plan_z) {
            const iy1 = (iy + offset)%yRange.size + y0;
            var elt = c_ptrTo(arr.localAccess[ix, iy1, z0]);
            plan_z.execute(elt, elt);
            ref srcRef = arr.localAccess[ix, iy1, z0];
            ref dstRef = dest[iy1,ix,z0];
            __primitive("chpl_comm_put", srcRef, dstRef.locale.id, dstRef, myLineSize);
          }
        }
        tt.stop(TimeStages.Execute);

        // on loc ends
      }
    }
  }

  // Set up the X FFT plan.
  // Assumes that we get the XZ plane.
  proc setupXPlan(param ftType : FFTtype, xzplane : [?Dom] ?T, signOrKind, flags : c_uint) {
    // Pull signOrKind locally since this may be an array
    // we need to take a pointer to.
    var mySignOrKind = signOrKind;
    var arg0 : _signOrKindType(ftType);
    select ftType {
        when FFTtype.R2R do arg0 = c_ptrTo(mySignOrKind);
        when FFTtype.DFT do arg0 = mySignOrKind;
      }

    // Set FFTW parameters
    const myDom = xzplane.localSubdomain();
    const xRange = Dom.dim(1);
    const zRange = myDom.dim(3);
    var nn = xRange.size : c_int;
    var nnp = c_ptrTo(nn);
    var howmany = myDom.dim(3).size : c_int;
    var rank = 1 : c_int;
    var stride = myDom.dim(3).size : c_int;
    var idist = 1 : c_int;
    var arr0 = c_ptrTo(xzplane.localAccess[myDom.first]);
    return new FFTWplan(ftType, rank, nnp, howmany, arr0,
                        nnp, stride, idist,
                        arr0, nnp, stride, idist,
                        arg0, flags);
  }

  /* X helper.

     See TODO for question on sign
  */
  proc doFFT_X(param ftType: FFTtype, arr : [?Dom] ?T, warmUpOnly : bool, signOrKind, flags : c_uint) {

    coforall loc in Locales {
      on loc {

        // Set up for the yz transforms on each domain.
        const myDom = arr.localSubdomain();
        const localIndex = myDom.first;
        const xRange = Dom.dim(1);
        const zRange = myDom.dim(3);

        // Split the y-range. Make this dimension agnostic
        const yChunk = chunk(myDom.dim(2), numLocales, here.id);
        const numberOfActualPlanes = if yChunk.size > numberOfPlanes then numberOfPlanes else yChunk.size;

        if (warmUpOnly) {
          var plan_x = setup1DPlan(T, ftType, xRange.size, zRange.size, signOrKind, flags);
        } else {
          // We assume numberOfActualPlanes tasks
          coforall iplane in 0.. #numberOfActualPlanes {
            var myplane : [{xRange, 0..0, zRange}] T;
            var plan_x = setup1DPlan(T, ftType, xRange.size, zRange.size, signOrKind, flags);
            var tt = new TimeTracker();

            var myChunk = chunk(yChunk, numberOfActualPlanes, iplane);
            const x0 = xRange.first;
            for j in myChunk {
              tt.start();
              myplane = arr[{xRange,j..j,zRange}];
              tt.stop(TimeStages.Comms);

              tt.start();
              forall iz in zRange with (ref plan_x) {
                var elt = c_ptrTo(myplane[x0, 0, iz]);
                plan_x.execute(elt, elt);
              }
              tt.stop(TimeStages.Execute);

              // Push back the pencil here
              tt.start();
              arr[{xRange,j..j,zRange}] = myplane;
              tt.stop(TimeStages.Comms);
            }
          }
        }
      }
      // End of on-loc
    }
  }

  /* X helper, transposed.

     FFT the x direction of arr, and store transposed in dest.

     x_y_z -> y_x_z
  */
  proc doFFT_X_Transposed(param ftType: FFTtype, 
                          dest : [?DomDest] ?T, warmUpOnly : bool, signOrKind, flags : c_uint) {

    coforall loc in Locales {
      on loc {

        // Split the y-range. Make this dimension agnostic
        // Transposed dimensions
        const myDomDest = dest.localSubdomain();
        const yRange = myDomDest.dim(1);
        const xRange = myDomDest.dim(2);
        const zRange = myDomDest.dim(3);

        const numTasks = min(here.maxTaskPar, zRange.size);
        const numTransforms = zRange.size/numTasks;
        var plan_x = setupPlanColumns(T, ftType, {xRange, zRange}, numTransforms, signOrKind, flags);
        var plan_x1 = setupPlanColumns(T, ftType, {xRange, zRange}, numTransforms+1, signOrKind, flags);

        if (!warmUpOnly) {
          const x0 = xRange.first;

          coforall itask in 0..#numTasks with (ref plan_x, ref plan_x1) {
            const myzRange = chunk(zRange, numTasks, itask);
            for iy in yRange {
              var elt = c_ptrTo(dest.localAccess[iy, x0, myzRange.first]);
              select myzRange.size {
                  when numTransforms do plan_x.execute(elt,elt);
                  when numTransforms+1 do plan_x1.execute(elt,elt);
                  otherwise halt("Bad myzRange size");
                }
            }
          }
        }
      }
    }
  }

  // I could not combine these, so keep them separate for now.
  private proc _signOrKindType(param ftType : FFTtype) type
    where (ftType==FFTtype.DFT) {
    return c_int;
  }
  private proc _signOrKindType(param ftType : FFTtype) type
    where (ftType==FFTtype.R2R) {
    return c_ptr(fftw_r2r_kind);
  }


  module FFT_Timers {
    use Time;
    // Time the various FFT steps.
    config const timeTrackFFT=false;

    enum TimeStages {X, YZ, Execute, Memcpy, Comms};
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
