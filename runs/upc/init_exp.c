#include <math.h>
#ifdef M_PI
#  define PI	M_PI
#elif !defined(PI)
#  define PI           3.14159265358979323846  
#endif

void init_exp(double *ex, double alpha, int THREADS, int MYTHREAD, int NX, int NY, int NZ, 
	      int TY, int TZ)
{
  double   kon;
  double   *exptr;
  int i,j,k,t,n;
  int it,jt,kt;
  int i_start, i_end, j_start, j_end;
  double ex_0 = 1.0;
  double ex_1;
  double prev;
  double ii,jj,kk;
  int myrow= MYTHREAD % TY;
  int myplane = MYTHREAD/TY;
  
  n=0;
  kon  = -4.*alpha*PI*PI;
  
  i_start = (NX/TY)*(myrow);
  i_end = (NX/TY)*(myrow+1);
  j_start = (NY/TZ)*(myplane);
  j_end = (NY/TZ)*(myplane+1);
	
  /*intialize the exp array*/
  /*by the time we use the exp array the layout will be 
    NY planes, NX rows, NZ columns
  */
  
  for(j=j_start; j < j_end; j++) {
    jt = (j >= NY/2 ? j-NY : j);
    jj = jt*jt;
    for(i=i_start; i < i_end; i++) {
      it = (i > NX/2 ? i-NX : i);
      ii = it*it;
      for (k=0; k < NZ; k++) {
        kt = (k >= NZ/2 ? k-NZ : k);
	kk = kt*kt;
	ex[n] = exp(kon*(ii + jj + kk));
	//	printf("%d> exp: %d %d %d %g\n", MYTHREAD, i, j, k, ex[n]);
	n++;
      }
    }
  }
}
