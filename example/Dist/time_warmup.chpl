use DistributedFFT;
use Time;

config const Ng=64;
config const ntrials=10;

const Dom = newSlabDom((64,64,64));
var arr : [Dom] complex;

var timeit : Timer;

timeit.clear(); timeit.start();
warmUpPlanner(arr);
timeit.stop();
writef("Time to warm up planner.... %r\n",timeit.elapsed());

var tmax = 0.0;
var tsum = 0.0;
for ii in 1..ntrials {
  timeit.clear(); timeit.start();
  warmUpPlanner(arr);
  timeit.stop();
  const t0 = timeit.elapsed();
  if (t0 > tmax) then tmax=t0;
  tsum += t0;
}


writef("Avg/Max time for planners : %r / %r\n",tsum/ntrials, tmax);

