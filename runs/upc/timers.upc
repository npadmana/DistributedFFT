#include <inttypes.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <timers.uph>
#include <upc.h>
#include <upc_collective.h>

static char *FFTimers_descr[T_NUMTIMERS] = {TIMER_STR_NAMES};
typedef uint64_t    ft_timer_t;

static ft_timer_t FTTimers_begin[T_NUMTIMERS];
static ft_timer_t FTTimers_total[T_NUMTIMERS] = { 0 };

static shared [T_NUMTIMERS] ft_timer_t alltimers[T_NUMTIMERS * THREADS];
static shared [] ft_timer_t *timers0;

/*
 * Timers
 */
#if defined(__UPC_TICK__)
#define get_ticks() upc_ticks_now()
#else
static inline
uint64_t get_ticks() {
  uint64_t retval;
  struct timeval tv;
  if (gettimeofday(&tv, NULL)) {
    perror("gettimeofday");
    abort();
  }
  retval = ((int64_t)tv.tv_sec) * 1000000 + tv.tv_usec;
  return retval;
}
#endif

void
timer_update(int tid, int action)
{
  if (action == FT_TIME_BEGIN) {
    FTTimers_begin[tid] = get_ticks();
  }
  else { /* FT_TIME_END */
    FTTimers_total[tid] += (get_ticks() - FTTimers_begin[tid]);
  }
}

void timer_clear()
{
  int i;
  for(i=0; i<T_NUMTIMERS; i++) {
    FTTimers_total[i] = 0;
  }
}

uint64_t timer_val(int tid)
{
#if defined(__UPC_TICK__)
  return (uint64_t) upc_ticks_to_ns(FTTimers_total[tid]) * 1000;
#else
  return FTTimers_total[tid];
#endif
}

char *timer_descr(int tid)
{
  return FFTimers_descr[tid];
}


void print_all_timers(char class, int TY, int TZ) {
  int timer=0, thread=0;
  
  upc_barrier;
  for(timer=0; timer<T_NUMTIMERS; timer++) {
    uint64_t timerval = timer_val(timer);
    alltimers[MYTHREAD * T_NUMTIMERS + timer] = timerval;
  }
  if (MYTHREAD == 0) {
    timers0 = (shared [] ft_timer_t *) upc_alloc(T_NUMTIMERS * sizeof(ft_timer_t)*THREADS);
  }
  upc_barrier;
  upc_all_gather(timers0, alltimers, sizeof(ft_timer_t)*T_NUMTIMERS, UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
  
  if(MYTHREAD==0) {
    for(thread=0; thread<THREADS; thread++) {
      printf("THREAD: %d COMM: UPC CLASS: %c THREADS: %d TY: %d TZ: %d ", thread, class, THREADS, TY, TZ);
      for(timer=0; timer<T_NUMTIMERS; timer++) {
	printf(" (%s) %g", FFTimers_descr[timer], ((double) timers0[timer+thread*T_NUMTIMERS])/1.e6);
      }
      printf("\n");
    }
    upc_free(timers0);
  }
  upc_barrier;
}
