#ifndef __MYFFT_H__
#define __MYFFT_H__ 1

#include <math.h>

typedef struct tagComplex_t {
	double real;
	double imag;
} myComplex_t;

typedef enum {MYFFT_FORWARD=1, MYFFT_BACKWARD=-1} myfftdir_t;

struct myfft_plan_t{
  int instride;
  int indist;
  int outstride;
  int outdist;
  myfftdir_t direction;
  int len;
  int numffts;
  int flags;
  void* lib_plan;
};

#define MYFFT_LIB_HANDLES_STRIDES 1<<1

void myfft_plan_1d(myfftdir_t direction, int len, int numffts,
		   myComplex_t *sample_outarray,
		   int outstride, int outdist,
		   myComplex_t *sample_inarray,
		   int instride, int indist,
		   struct myfft_plan_t *plan,
		   int flags);

void run_1d_fft(myComplex_t * out, myComplex_t *in, struct myfft_plan_t *plan);

char *myfft_id_str();

#endif
