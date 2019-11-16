#include <fftw3.h>
#include <stdlib.h>
#include <myfft.h>
#include <stdio.h>
#include <string.h>

char *myfft_id_str(char *buffer) {
  strcpy(buffer, "FFTW3");
  return buffer;
}

typedef struct tagFFTW_PLAN{
  fftw_plan *plan;
  myComplex_t *buffer;
}myfftw_plan_t;

void myfft_plan_1d(myfftdir_t direction, int len, int numffts, 
		   myComplex_t *sample_outarray,
		   int outstride, int outdist, 
		   myComplex_t *sample_inarray,
		   int instride, int indist,
		   struct myfft_plan_t *plan, int flags) {
  
  myfftw_plan_t *myplan;

  myplan = malloc(sizeof(myfftw_plan_t));
  myplan->plan = malloc(sizeof(fftw_plan));
  
  if(flags & MYFFT_LIB_HANDLES_STRIDES) {
    *(myplan->plan) = 
      fftw_plan_many_dft(1, &len, numffts,
			 (fftw_complex*) sample_inarray, NULL, instride, indist,
			 (fftw_complex*) sample_outarray, NULL, outstride, outdist, 
			 (direction == MYFFT_FORWARD ? FFTW_FORWARD : FFTW_BACKWARD),
			 FFTW_MEASURE);
  } else {
    myplan->buffer = fftw_malloc(sizeof(myComplex_t)*len);
    *(myplan->plan) = 
      fftw_plan_many_dft(1, &len, 1,
			 (fftw_complex*) myplan->buffer, NULL, 1, len,
			 (fftw_complex*) myplan->buffer, NULL, 1, len, 
			 (direction == MYFFT_FORWARD ? FFTW_FORWARD : FFTW_BACKWARD),
			 FFTW_MEASURE);
  }
  
  plan->instride = instride;
  plan->outstride =outstride;
  plan->indist = indist;
  plan->outdist = outdist;
  plan->direction = direction;
  plan->len = len;
  plan->numffts = numffts;
  plan->lib_plan = myplan;
  plan->flags = flags;  
}

void run_1d_fft(myComplex_t * out, myComplex_t *in, struct myfft_plan_t *plan) {
  fftw_plan *lib_plan = (fftw_plan*) ((myfftw_plan_t*)(plan->lib_plan))->plan;
  
  if(plan->flags & MYFFT_LIB_HANDLES_STRIDES) {
    fftw_execute_dft(*lib_plan, (fftw_complex *)in, (fftw_complex *)out);
  } else {
    myComplex_t *buffer = ((myfftw_plan_t*) (plan->lib_plan))->buffer;
    int i,j;
    myComplex_t *temp_in=in, *temp_out=out;
    for(i=0; i<plan->numffts; i++) {
      for(j=0; j<plan->len; j++) {
	buffer[j] = in[j*plan->instride];
      }
      fftw_execute_dft(*lib_plan, (fftw_complex *)buffer, (fftw_complex *)buffer);
      for(j=0; j<plan->len; j++) {
	out[j*plan->outstride] = buffer[j];
      }
      in+=plan->indist;
      out+=plan->outdist;
    }
  }
}

