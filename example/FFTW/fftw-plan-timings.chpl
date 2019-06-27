use Time;
use FFTW;

var timeit : Timer;


var fftw_planner_lock$ : sync bool;


config const Ng=256;
config const ntry=1000;

const Dom = {0.. #Ng, 0.. #Ng, 0.. #Ng};
var arr : [Dom] complex;

// Warmup plans
timeit.clear(); timeit.start();
var plan = plan_dft(arr, FFTW_FORWARD, FFTW_MEASURE);
timeit.stop();
writef("Elapsed time for plan warmup : %r \n", timeit.elapsed());

// Now do this in a forall loop
timeit.clear(); timeit.start();
var sum : int = 0;
forall ii in 1..ntry with (+ reduce sum,
                           var locarr : [Dom] complex,
                           var myplan = new fftw_plan_rec(locarr, FFTW_FORWARD, FFTW_MEASURE)) {
  sum += 1;
}
timeit.stop();
writef("Sum= %i\n",sum);
writef("Max time for plan creation in loop: %r \n", timeit.elapsed());

record fftw_plan_rec {
  var plan : fftw_plan;

  proc init(arr, sign, flags) {
    fftw_planner_lock$.writeEF(true);
    plan = plan_dft(arr, sign, flags);
    fftw_planner_lock$.readFE();
  }

  proc deinit() {
    destroy_plan(plan);
  }

  proc execute() {
    execute(plan);
  }
}

