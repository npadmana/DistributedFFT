/* Documentation for DistributedFFT */
prototype module DistributedFFT {

  use BlockDist;
  use AllLocalesBarriers;
  use RangeChunk;
  use FFTW;
  use FFTW.C_FFTW;
  use FFT_Locks;
  use FFT_Timers;
  require "npFFTW.h";

  config param usePerformant=true;
  config param usePrimitiveComm=true;

  proc deinit() {
    cleanup();
  }

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
  // https://github.com/chapel-lang/chapel/issues/13319
  pragma "default intent is ref"
  record FFTWplan {
    param ftType : FFTtype;
    var plan : fftw_plan;

    // Mimic the advanced interface 
    proc init(param ftType : FFTtype, args ...?k) {
      this.ftType = ftType;
      this.complete();
      plannerLock.lock();
      select ftType {
        when FFTtype.DFT do plan = fftw_plan_many_dft((...args));
        when FFTtype.R2R do plan = fftw_plan_many_r2r((...args));
      }
      plannerLock.unlock();
    }

    proc deinit() {
      plannerLock.lock();
      destroy_plan(plan);
      plannerLock.unlock();
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
      execute(c_ptrTo(arr1), c_ptrTo(arr2));
    }

    inline proc execute(ref arr1 : ?T) where (!isAnyCPtr(T)) {
      execute(arr1, arr1);
    }

    proc isValid : bool {
      extern proc isNullPlan(plan : fftw_plan) : c_int;
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
    if (usePerformant) {
      doFFT_Transposed_Performant(ftType, src, dest, signOrKind);
    } else {
      doFFT_Transposed_Naive(ftType, src, dest, signOrKind);
    }
  }


  /* FFT.

     Stores the FFT in Dst transposed (xyz -> yxz).
   */
  proc doFFT_Transposed_Naive(param ftType : FFTtype,
                              Src: [?SrcDom] ?T,
                              Dst : [?DstDom] T,
                              signOrKind) {
    // Sanity checks
    if SrcDom.rank != 3 || DstDom.rank != 3 then compilerError("Code is designed for 3D arrays only");
    if SrcDom.dim(1) != DstDom.dim(2) then halt("Mismatched x-y ranges");
    if SrcDom.dim(2) != DstDom.dim(1) then halt("Mismatched y-x ranges");
    if SrcDom.dim(3) != DstDom.dim(3) then halt("Mismatched z ranges");

    coforall loc in Locales do on loc {
      var timeTrack = new TimeTracker();

      const (xSrc, ySrc, zSrc) = SrcDom.localSubdomain().dims();
      const (yDst, xDst, _) = DstDom.localSubdomain().dims();

      // Set up FFTW plans
      var xPlan = setup1DPlan(T, ftType, xDst.size, zSrc.size, signOrKind, FFTW_MEASURE);
      var yPlan = setup1DPlan(T, ftType, ySrc.size, zSrc.size, signOrKind, FFTW_MEASURE);
      var zPlan = setup1DPlan(T, ftType, zSrc.size, 1, signOrKind, FFTW_MEASURE);

      // Use temp work array to avoid overwriting the Src array
      var myplane : [{0..0, ySrc, zSrc}] T;

      for ix in xSrc {
        // Copy source to temp array
        myplane = Src[{ix..ix, ySrc, zSrc}];

        // Y-transform
        timeTrack.start();
        forall iz in zSrc {
          yPlan.execute(myplane[0, ySrc.first, iz]);
        }
        timeTrack.stop(TimeStages.Y);

        // Z-transform, offset to reduce comm congestion/collision
        timeTrack.start();
        forall iy in offset(ySrc) {
          zPlan.execute(myplane[0, iy, zSrc.first]);
          // Transpose data into Dst
          Dst[{iy..iy, ix..ix, zSrc}] = myplane[{0..0, iy..iy, zSrc}];
        }
        timeTrack.stop(TimeStages.Z);
      }

      // Wait until all communication is complete
      allLocalesBarrier.barrier();

      // X-transform
      timeTrack.start();
      forall (iy, iz) in {yDst, zSrc} {
        xPlan.execute(Dst[iy, xDst.first, iz]);
      }
      timeTrack.stop(TimeStages.X);
    }
  }


  /* FFT.

     Stores the FFT in Dst transposed (xyz -> yxz).
   */
  proc doFFT_Transposed_Performant(param ftType : FFTtype,
                                   Src: [?SrcDom] ?T,
                                   Dst : [?DstDom] T,
                                   signOrKind) {
    // Sanity checks
    if SrcDom.rank != 3 || DstDom.rank != 3 then compilerError("Code is designed for 3D arrays only");
    if SrcDom.dim(1) != DstDom.dim(2) then halt("Mismatched x-y ranges");
    if SrcDom.dim(2) != DstDom.dim(1) then halt("Mismatched y-x ranges");
    if SrcDom.dim(3) != DstDom.dim(3) then halt("Mismatched z ranges");

    coforall loc in Locales do on loc {
      var timeTrack = new TimeTracker();

      const (xSrc, ySrc, zSrc) = SrcDom.localSubdomain().dims();
      const (yDst, xDst, _) = DstDom.localSubdomain().dims();
      const myLineSize = zSrc.size*numBytes(T);

      // Setup FFTW plans
      var yPlan = setupBatchPlanColumns(T, ftType, {ySrc, zSrc}, parDim=2, signOrKind, FFTW_MEASURE);
      var xPlan = setupBatchPlanColumns(T, ftType, {xDst, zSrc}, parDim=2, signOrKind, FFTW_MEASURE);
      var zPlan = setup1DPlan(T, ftType, zSrc.size, 1, signOrKind, FFTW_MEASURE);

      // Use temp work array to avoid overwriting the Src array
      var myplane : [{0..0, ySrc, zSrc}] T;

      if usePrimitiveComm {
        forall iy in ySrc {
          copy(myplane[0, iy, zSrc.first], Src[xSrc.first, iy, zSrc.first], myLineSize);
        }
      } else {
        myplane = Src[{xSrc.first..xSrc.first, ySrc, zSrc}];
      }

      for ix in xSrc {
        // Y-transform
        timeTrack.start();
        forall (plan, myzRange) in yPlan.batch() {
          plan.execute(myplane[0, ySrc.first, myzRange.first]);
        }
        timeTrack.stop(TimeStages.Y);

        // Z-transform, offset to reduce comm congestion/collision
        timeTrack.start();
        forall iy in offset(ySrc) {
          zPlan.execute(myplane[0, iy, zSrc.first]);
          // Transpose data into Dst, and copy the next Src slice into myplane
          if usePrimitiveComm {
            copy(Dst[iy, ix, zSrc.first], myplane[0, iy, zSrc.first], myLineSize);
          } else {
            Dst[{iy..iy, ix..ix, zSrc}] = myplane[{0..0, iy..iy, zSrc}];
          }
          if (ix != xSrc.last) {
            if usePrimitiveComm {
              copy(myplane[0, iy, zSrc.first], Src[ix+1, iy, zSrc.first], myLineSize);
            } else {
              myplane[{0..0, iy..iy, zSrc}] = Src[{ix+1..ix+1, iy..iy, zSrc}];
            }
          }
        }
        timeTrack.stop(TimeStages.Z);
      }

      // Wait until all communication is complete
      allLocalesBarrier.barrier();

      // X-transform
      timeTrack.start();
      forall (plan, myzRange) in xPlan.batch() {
        for iy in yDst {
          plan.execute(Dst[iy, xDst.first, myzRange.first]);
        }
      }
      timeTrack.stop(TimeStages.X);
    }
  }

  iter offset(r: range) { halt("Serial offset not implemented"); }
  iter offset(param tag: iterKind, r: range) where (tag==iterKind.standalone) {
    forall i in r + (r.size/numLocales * here.id) do {
      yield i % r.size + r.first;
    }
  }

  proc copy(ref dst, const ref src, numBytes: int) {
    if dst.locale.id == here.id {
      __primitive("chpl_comm_get", dst, src.locale.id, src, numBytes.safeCast(size_t));
    } else if src.locale.id == here.id {
      __primitive("chpl_comm_put", src, dst.locale.id, dst, numBytes.safeCast(size_t));
    } else {
      halt("Remote src and remote dst not yet supported");
    }
  }


  pragma "default intent is ref"
  record BatchedFFTWplan {
    param ftType : FFTtype;
    const parRange: range;
    const numTasks: int;
    const batchSizeSm, batchSizeLg: int;
    var planSm, planLg: FFTWplan(ftType);

    proc init(type arrType, param ftType : FFTtype, dom : domain(2), parDim : int, signOrKind, in flags : c_uint) {
      this.ftType = ftType;
      this.parRange = dom.dim(parDim);
      this.numTasks = min(here.maxTaskPar, parRange.size);
      this.batchSizeSm = parRange.size/numTasks;
      this.batchSizeLg = parRange.size/numTasks+1;
      this.planSm = setupPlanColumns(arrType, ftType, dom, batchSizeSm, signOrKind, flags);
      this.planLg = setupPlanColumns(arrType, ftType, dom, batchSizeLg, signOrKind, flags);
    }

    iter batch() {
      halt("Serial iterator not implemented");
    }

    iter batch(param tag : iterKind) where (tag==iterKind.standalone) {
      coforall chunk in chunks(parRange, numTasks) {
        if chunk.size == batchSizeSm then yield (planSm, chunk);
        if chunk.size == batchSizeLg then yield (planLg, chunk);
      }
    }
  }

  proc setupBatchPlanColumns(type arrType, param ftType : FFTtype, dom : domain(2), parDim : int, signOrKind, in flags : c_uint) {
    return new BatchedFFTWplan(arrType, ftType, dom, parDim, signOrKind, flags);
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

  module FFT_Locks {
    // https://github.com/chapel-lang/chapel/issues/9881
    // https://github.com/chapel-lang/chapel/issues/12300
    private use ChapelLocks;
    pragma "locale private"
    var plannerLock : chpl_LocalSpinlock;
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
