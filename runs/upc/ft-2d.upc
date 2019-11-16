
#include <upc.h>
#include <fft3d.uph>
#include <stdio.h>
#include <params.h>
#include <stdlib.h>

#include "timers.uph"

void init_seed(double seed);
double randlc (double *x, double a);
void vranlc(int n, double *x, double a, double *y);
double ipow46(double, int, int);
void init_exp(double *ex, double alpha, int total, int who, int NX, int NY, int NZ, int TY, int TZ);
void checksum_verify (int d1, int d2, int d3, int nt, myComplex_t *sums);

shared myComplex_t *arr_fwd, *arr_back, *arrB, *arrC;
myComplex_t *myarr_fwd, *myarr_back, *myB, *myC;
double * exp_array;

unsigned int NX, NY, NZ;
int TY, TZ;
unsigned int NX_TY, NX_TZ;
unsigned int NY_TY, NY_TZ;
unsigned int NZ_TY, NZ_TZ;

int myplane, myrow;
unsigned int NX_start, NX_end;
unsigned int NY_start, NY_end;
unsigned int NZ_start, NZ_end;

unsigned int max_size;

#define PAD(v) (v)

void parabolic2(myComplex_t *dest, myComplex_t *src) {
  int i,j,k,idx;
  TIMER_START(T_EVOLVE);
  for(k=0; k<NY_TZ; k++) {
    for(j=0; j<NX_TY; j++) {
      for(i=0; i<NZ; i++) {
	idx = k*PAD(NX_TY)*PAD(NZ)+j*PAD(NZ)+i;
	//dest[idx].real =src[idx].real+1.0;
	//dest[idx].imag = src[idx].imag+1.0;
	dest[idx].real = src[idx].real=src[idx].real*exp_array[k*(NX_TY)*NZ+j*NZ+i];
	dest[idx].imag = src[idx].imag=src[idx].imag*exp_array[k*(NX_TY)*NZ+j*NZ+i];
      }
    }
  }
  TIMER_STOP(T_EVOLVE);
}

shared myComplex_t allsums[THREADS];
myComplex_t checksum(myComplex_t *arr) {
  myComplex_t sum, ret;
  int q,r,s,proc,j;
  double dNX, dNY, dNZ;
  double dNTOTAL;	
  TIMER_START(T_CHECKSUM);
  dNX = (double) NX;
  dNY = (double) NY;
  dNZ = (double) NZ;
  dNTOTAL = (dNX*dNY*dNZ);
  sum.real = 0;
  sum.imag = 0;
  for(j=1; j<=1024; j++) {
    q = j % NX;
    r = (3*j) % NY;
    s = (5*j) % NZ;
    
    /*at this point planes are oriented
      NX planes, NZ rows, NY columns*/
    proc = (q/(NX_TZ))*TY +  (s/(NZ_TY));
    
    if(MYTHREAD==proc) {
      myComplex_t elem = arr[(q%(NX_TZ))*(PAD(NY)*PAD(NZ_TY)) + r + (s%(NZ_TY))*PAD(NY)];
      sum.real+=elem.real;
      sum.imag+=elem.imag;
    }
  }
  
  allsums[MYTHREAD] = sum;
  upc_barrier;
  if(MYTHREAD==0) {
    double allreal=0, allimag=0;
    int t=0;
    for(t=0; t<THREADS; t++) {
      allreal+=allsums[t].real;
      allimag+=allsums[t].imag;
    }
    ret.real = allreal;
    ret.imag = allimag;
  }

  if(MYTHREAD==0) {
    ret.real /=dNTOTAL;
    ret.imag /=dNTOTAL;
  }
  TIMER_STOP(T_CHECKSUM);
  TIMER_START(T_BARRIER_CHK);
  upc_barrier;
  TIMER_STOP(T_BARRIER_CHK);
  return ret;
}

void allocate_arrays() {

  max_size = NX*(NY_TY)*(NZ_TZ);

  arr_fwd = (shared myComplex_t*) upc_all_alloc(THREADS, sizeof(myComplex_t)*max_size);
  arr_back = (shared myComplex_t*) upc_all_alloc(THREADS, sizeof(myComplex_t)*max_size);
  arrB = (shared myComplex_t*) upc_all_alloc(THREADS, sizeof(myComplex_t)*max_size);
  arrC = (shared myComplex_t*) upc_all_alloc(THREADS, sizeof(myComplex_t)*max_size);
  
  exp_array = (double*) malloc(sizeof(double)*(NX_TY)*(NY_TZ)*NZ);
  myarr_fwd = (myComplex_t*) ((shared [] myComplex_t*) &arr_fwd[MYTHREAD]);
  myarr_back = (myComplex_t*) ((shared [] myComplex_t*) &arr_back[MYTHREAD]);
  myB = (myComplex_t*) ((shared [] myComplex_t*) &arrB[MYTHREAD]);
  myC = (myComplex_t*) ((shared [] myComplex_t*) &arrC[MYTHREAD]);
  
  myplane = MYTHREAD / TY;
  myrow = MYTHREAD % TY;
  
  NX_start = 0; NX_end = NX-1;
  NY_start = myrow*(NY_TY); NY_end = (myrow+1)*(NY_TY)-1;
  NZ_start = myplane*(NZ_TZ); NZ_end = (myplane+1)*(NZ_TZ)-1;
}


void init_arrays(myComplex_t *arr, myComplex_t *scratch) {

  /*fill in later with real values*/
  /*for now just initialze the array w/ its index*/
  int i, j, k;
  double  x0, start, an;
	
  TIMER_START(T_SETUP);
  start = SEED;
  
  init_seed(start);
  an  = ipow46(A, 2*NX, NZ_start*NY+NY_start);
  randlc(&start, an);
  an = ipow46(A, 2*NX, NY);
  
  for(k=0; k<NZ_TZ; k++) {
    x0 = start;
    vranlc(2*NX*(NY_TY), &x0, A, (double*) (arr+k*(NX)*(NY_TY)));
    if(k!=(NZ_TZ-1)) randlc(&start,an);
  }

  init_exp(exp_array, 1.0e-6, THREADS, MYTHREAD, NX, NY, NZ, TY, TZ);	
  TIMER_STOP(T_SETUP);
  
  TIMER_START(T_BARRIER_WAIT);
  upc_barrier;
  TIMER_STOP(T_BARRIER_WAIT);
}

void dumpArray(char *str, myComplex_t *in) {
  int i,j,k;
  upc_barrier;
  
  for(j=0; j<max_size; j++) {
    fprintf(stdout, "%d> %s (idx: %d) %g %g\n", MYTHREAD, str, j, in[j].real, in[j].imag);
  }
  upc_barrier;
}

void print_usage(char *str) {
  upc_barrier;
  if(MYTHREAD==0) fprintf(stderr, "useage: %s TY TZ\n", str);
  upc_barrier;
}

int main(int argc, char **argv) {
  
  int blocking_factor_x, bfx_set=0;
  int blocking_factor_y, bfy_set=0;
  int iter;
  int i, j;
  int PX, PY, PZ, PT;
  myComplex_t checksums[MAX_ITER];
  fft3d_handle_t *fft_handle_forward, *fft_handle_backward;
  
  NX=NX_IN;
  NY=NY_IN;
  NZ=NZ_IN;
  
  if(argc<3) {print_usage(argv[0]); exit(1);};

  PX = atoi(argv[1]);
  PY = 1 ;
  PZ = atoi(argv[2]);
  PT = 1;
  
  TY = PX*PY;
  TZ = PZ*PT;
  
  NX_TY = NX/TY; NX_TZ = NX/TZ;
  NY_TY = NY/TY; NY_TZ = NY/TZ;
  NZ_TY = NZ/TY; NZ_TZ = NZ/TZ;
  
  if(TY*TZ!=THREADS) { 
    if(MYTHREAD==0) fprintf(stderr, "Bad Arguments: TY*TZ!=THREADS\n"); exit(1);
  }
  if(TZ>NZ) { 
    if(MYTHREAD==0) fprintf(stderr, "Bad Arguments: TZ>NZ\n"); exit(1);
  }
  
  if(MYTHREAD==0) {
    fprintf(stderr, "starting %s (NX:%d,NY:%d,NZ:%d) on %d x %d processor grid\n", class_id_str, NX, NY, NZ, TY, TZ);
  }

  allocate_arrays();
  timer_clear();
  init_arrays(myarr_fwd, myB);
  
  fft_handle_forward  = plan3DFFT(NX, NY, NZ, PX, PY, PZ, PT, NY_TY, NX_TY, myB, myC, MYFFT_FORWARD);
  fft_handle_backward = plan3DFFT(NZ, NX, NY, PX, PY, PZ, PT, NX_TY, NZ_TY, myB, myC, MYFFT_BACKWARD);
  
  upc_barrier;
  
  run3DFFT(fft_handle_forward, arrB, arr_fwd, arrC);
  timer_clear();

  TIMER_START(T_TOTAL);
  init_arrays(myarr_fwd, myB);
  
  upc_barrier;
  run3DFFT(fft_handle_forward, arrB, arr_fwd, arrC);
  
  for(iter=1; iter <= MAX_ITER; iter++) {
    parabolic2(myarr_fwd, myB);
    run3DFFT(fft_handle_backward, arr_back, arr_fwd, arrC);
    checksums[iter-1] = checksum(myarr_back);
    
    if (MYTHREAD==0)
      fprintf(stdout, " 0> %30s %2d: %#17.14g %#17.14g\n",
	      "Checksum", iter, checksums[iter-1].real, checksums[iter-1].imag);
    fflush(stdout);
  }
  TIMER_STOP(T_TOTAL);
  
  if(MYTHREAD==0) {
    checksum_verify(NX,NY,NZ, MAX_ITER, checksums);
  }
  upc_barrier;
  
  if(MYTHREAD==0){
    double mflops;
    double mflop_rate;
    char idbuffer[10];
    char idbuffer2[20];
    
    mflops = (double)(1.0e-6*((double)NX*NY*NZ) *
		      (14.8157 + 7.19641*log((double)NX*NY*NZ)
		       + (5.23518 + 7.21113*log((double)NX*NY*NZ))*MAX_ITER));
    mflop_rate = mflops/(((double)timer_val(T_TOTAL))/1e6);
    printf("Total running time is %f s\n", ((double)timer_val(T_TOTAL))/1e6);
    printf("NAS FT (BACKEND: %s) (UPC: %s) CLASS (%c) on %d processors TY: %d TZ: %d Mflops: %g MFlop/s MFlop/s/Thread: %g\n", 
	    myfft_id_str(idbuffer), fft3d_id_str(idbuffer2), class_id_char, 
	    THREADS, TY, TZ, mflop_rate, mflop_rate/THREADS);
  }
  upc_barrier;
/*  
  print_all_timers(class_id_char, TY, TZ);
*/  
  upc_barrier;
  if(MYTHREAD==0){ fprintf(stderr, "benchmark finished calling final barrier\n");}
  upc_barrier;
  
  return 0;
}
