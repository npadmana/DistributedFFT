/* Documentation for DistributedFFT */
prototype module DistributedFFT {

  use BlockDist;
  use AllLocalesBarriers;
  use ReplicatedVar;
  use RangeChunk;
  use FFTW;
  use FFTW.C_FFTW;
  use FFT_Timers;
  require "npFFTW.h";

  config const useElegant=false;

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
      FFTW.execute(plan);
    }

    proc execute(arr1 : c_ptr(?T), arr2 : c_ptr(T)) {
      select ftType {
          when FFTtype.DFT do fftw_execute_dft(plan, arr1, arr2);
          when FFTtype.R2R do fftw_execute_r2r(plan, arr1, arr2);
        }
    }

    inline proc execute(ref arr1 : ?T, ref arr2 : T) where (!isAnyCPtr(T)) {
      var elt1 = c_ptrTo(arr1);
      var elt2 = c_ptrTo(arr2);
      select ftType {
          when FFTtype.DFT do fftw_execute_dft(plan, elt1, elt2);
          when FFTtype.R2R do fftw_execute_r2r(plan, elt1, elt2);
        }
    }

    inline proc execute(ref arr1 : ?T) where (!isAnyCPtr(T)) {
      var elt1 = c_ptrTo(arr1);
      select ftType {
          when FFTtype.DFT do fftw_execute_dft(plan, elt1, elt1);
          when FFTtype.R2R do fftw_execute_r2r(plan, elt1, elt1);
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

  proc doFFT_Transposed(param ftType : FFTtype,
                        src: [?SrcDom] ?T,
                        dest : [?DestDom] T,
                        signOrKind) {
    if (useElegant) {
      doFFT_Transposed_Elegant(ftType, src, dest, signOrKind);
    } else {
      doFFT_Transposed_Performant(ftType, src, dest, signOrKind);
    }
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

        // Start tracking time
        var tt = new TimeTracker();

        for ix in xSrc {
          // Copy data over. This is a completely
          // local operation, so use localAccess.
          // Or a series of memcpys for performance.
          [(tmp, iy,iz) in myplane.domain] myplane[0, iy,iz] = src.localAccess[ix,iy,iz];

          // y-transform
          tt.start();
          const y0 = ySrc.first;
          forall iz in zSrc with (ref plan_y) {
            plan_y.execute(myplane[0,y0,iz]);
          }
          tt.stop(TimeStages.Y);

          // z-transform
          tt.start();
          const z0 = zSrc.first;
          // Offset to reduce collisions
          const offset = (ySrc.size/numLocales)*here.id;
          forall iy in 0.. #ySrc.size with (ref plan_z) {
            const iy1 = (iy + offset)%ySrc.size + y0;
            plan_z.execute(myplane[0,iy1,z0]);
            // This is the transpose step
            dest[{iy1..iy1,ix..ix,zSrc}] = myplane[{0..0, iy1..iy1,zSrc}];
          }
          tt.stop(TimeStages.Z);
        }

        // Wait until all communication is complete
        allLocalesBarrier.barrier();

        // x-transform
        tt.start();
        const x0 = xDest.first;
        forall (iy,iz) in {yDest, zSrc} with (ref plan_x) {
          plan_x.execute(dest.localAccess[iy, x0, iz]);
        }
        tt.stop(TimeStages.X);

        // End of on-loc
      }
    }


    // End of doFFT
  }

  /* FFT.

     Stores the FFT in dest transposed (xyz -> yxz).
   */
  proc doFFT_Transposed_Performant(param ftType : FFTtype,
                                   src: [?SrcDom] ?T,
                                   dest : [?DestDom] T,
                                   signOrKind) {
    // Sanity checks
    if SrcDom.rank != 3 then halt("Code is designed for 3D arrays only");
    if DestDom.rank != 3 then halt("Code is designed for 3D arrays only");
    if SrcDom.dim(1) != DestDom.dim(2) then halt("Mismatched x-y ranges");
    if SrcDom.dim(2) != DestDom.dim(1) then halt("Mismatched y-x ranges");
    if SrcDom.dim(3) != DestDom.dim(3) then halt("Mismatched z ranges");


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

        const numTasks = min(here.maxTaskPar, zSrc.size);
        const numTransforms = zSrc.size/numTasks;

        // Set up FFTW plans
        var plan_y = setupPlanColumns(T, ftType, {ySrc, zSrc}, numTransforms, signOrKind, FFTW_MEASURE);
        var plan_y1 = setupPlanColumns(T, ftType, {ySrc, zSrc}, numTransforms+1, signOrKind, FFTW_MEASURE);
        var plan_x = setupPlanColumns(T, ftType, {xDest, zSrc}, numTransforms, signOrKind, FFTW_MEASURE);
        var plan_x1 = setupPlanColumns(T, ftType, {xDest, zSrc}, numTransforms+1, signOrKind, FFTW_MEASURE);
        var plan_z = setup1DPlan(T, ftType, zSrc.size, 1, signOrKind, FFTW_MEASURE);

        // Use this as a temporary work array to
        // avoid rewriting the src array
        var myplane : [{0..0, ySrc, zSrc}] T;
        const z0 = zSrc.first;
        forall iy in ySrc {
          c_memcpy(c_ptrTo(myplane[0,iy,z0]),
                   c_ptrTo(src.localAccess[xSrc.first,iy,z0]),
                   myLineSize);
        }

        // Time tracking
        var tt = new TimeTracker();

        for ix in xSrc {
          // y-transform
          tt.start();
          const y0 = ySrc.first;
          forall myzRange in batchedRange(zSrc) with (ref plan_y, ref plan_y1) {
            select myzRange.size {
                when numTransforms do plan_y.execute(myplane[0,y0,myzRange.first]);
                when numTransforms+1 do plan_y1.execute(myplane[0,y0,myzRange.first]);
                otherwise halt("Bad myzRange size");
              }
          }
          tt.stop(TimeStages.Y);

          // z-transform
          // Offset to reduce collisions
          tt.start();
          const offset = (ySrc.size/numLocales)*here.id;
          forall iy in 0.. #ySrc.size with (ref plan_z) {
            const iy1 = (iy + offset)%ySrc.size + y0;
            plan_z.execute(myplane[0,iy1,z0]);
            // This is the transpose step
            remotePut(dest[iy1,ix,z0], myplane[0,iy1,z0], myLineSize);
            // If not last slice, copy over
            if (ix != xSrc.last) {
              c_memcpy(c_ptrTo(myplane[0,iy1,z0]),
                       c_ptrTo(src.localAccess[ix+1,iy1,z0]),
                       myLineSize);
            }
          }
          tt.stop(TimeStages.Z);
        }

        // Wait until all communication is complete
        allLocalesBarrier.barrier();

        // x-transform
        tt.start();
        const x0 = xDest.first;
        forall myzRange in batchedRange(zSrc) with (ref plan_x, ref plan_x1) {
          for iy in yDest {
            ref elt = dest.localAccess[iy, x0, myzRange.first];
            select myzRange.size {
                when numTransforms do plan_x.execute(elt);
                when numTransforms+1 do plan_x1.execute(elt);
                otherwise halt("Bad myzRange size");
              }
          }
        }
        tt.stop(TimeStages.X);

        // End of on-loc
      }
    }


    // End of doFFT
  }

  inline proc remotePut(ref dstRef, ref srcRef, numBytes : int) {
    __primitive("chpl_comm_put", srcRef, dstRef.locale.id, dstRef, numBytes);
  }

  iter batchedRange(r : range) {
    halt("Serial iterator not implemented");
  }

  iter batchedRange(param tag : iterKind, r : range)
    where (tag==iterKind.standalone)
  {
    const numTasks = min(here.maxTaskPar, r.size);
    coforall itask in 0.. #numTasks {
      const myr = chunk(r, numTasks, itask);
      yield myr;
    }
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

    enum TimeStages {X, Y, Z};
    const stageDomain = {TimeStages.X..TimeStages.Z};
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
