#include "fft3d.uph"
#include "timers.uph"
#include <stdlib.h>
#include <string.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define PAD(VAR) ((VAR)+0)

#ifndef CYCLE_EXCHANGE
#define CYCLE_EXCHANGE 0
#endif

#define MY_UPC_POLL() do{TIMER_START(T_POLL); TIMER_STOP(T_POLL);} while(0)
#define MY_UPC_BARRIER() do{TIMER_START(T_BARRIER_WAIT); upc_barrier; TIMER_STOP(T_BARRIER_WAIT);} while(0)

char *fft3d_id_str(char *buf) {
  sprintf(buf, "slabs");
  return buf;
}

#define LOCATION_TO_THREAD(handle, inX,inY,inZ,inT) ((inT)*(handle->PX)*(handle->PY)*(handle->PZ) + (inZ)*(handle->PY)*(handle->PX) + (inY)*(handle->PX) + inX)

fft3d_handle_t *plan3DFFT(unsigned int NX, unsigned int NY, unsigned int NZ, int PX, int PY, int PZ, int PT,
			  unsigned int blocking_factor_x, 
			  unsigned int blocking_factor_y, 
			  myComplex_t *sample_arrayA, 
			  myComplex_t* sample_arrayB, myfftdir_t fftdir) {
  fft3d_handle_t *ret;
	
  ret = (fft3d_handle_t*) malloc(sizeof(fft3d_handle_t));
  
  ret->NX = NX;
  ret->NY = NY;
  ret->NZ = NZ;
  
  ret->TY = PX*PY;
  ret->TZ = PZ*PT;
  
  ret->PX = PX;
  ret->PY = PY;
  ret->PZ = PZ;
  ret->PT = PT;
  
  ret->myloc[0] = MYTHREAD % PX;
  ret->myloc[1] = (MYTHREAD/PX) % PY;
  ret->myloc[2] = (MYTHREAD/(PX*PY)) % PZ;
  ret->myloc[3] = (MYTHREAD/(PX*PY*PZ)) % PT;
  
  ret->NTOTAL = NX*NY*NZ;
  ret->mythreadplane = MYTHREAD / ret->TY;
  ret->mythreadrow = MYTHREAD % ret->TY;
  
  ret->num_local_planes = NZ / ret->TZ;
  ret->num_local_rows = NY / ret->TY;
  
  blocking_factor_x = NY / ret->TY;
  blocking_factor_y = NX / ret->TY;

  ret->blocking_factor_x = blocking_factor_x;
  ret->blocking_factor_y = blocking_factor_y;
  
  {
    /*take max of number of msgs in round 1 v. round 2*/
    unsigned int max_num_msgs = MAX((NZ/ret->TZ)*((NY/ret->TY)/blocking_factor_x)*ret->TY /*# in round 1*/,
			      (NZ/ret->TZ)*((NX/ret->TY)/blocking_factor_y)*ret->TZ /*# in round 2*/);
    
    ret->comm_handles = malloc(sizeof(upc_handle_t)*max_num_msgs);
    ret->comm_handle_idx = 0;
  }
  
  myfft_plan_1d(fftdir, NX, blocking_factor_x, sample_arrayB, blocking_factor_x, 1, sample_arrayA,
                1, PAD(NX), &ret->plans[0], MYFFT_LIB_HANDLES_STRIDES);
  myfft_plan_1d(fftdir, NY, blocking_factor_y, sample_arrayB, blocking_factor_y, 1, sample_arrayA,
                1, PAD(NY), &ret->plans[1], MYFFT_LIB_HANDLES_STRIDES);
  myfft_plan_1d(fftdir, NZ, (NX/ret->TY), sample_arrayB, 1, PAD(NZ), sample_arrayA, 
		(NY/ret->TZ)*(NX/ret->TY), 1,
                &ret->plans[2], MYFFT_LIB_HANDLES_STRIDES);
  return ret;
}


#define GET_LOCAL_PORTION(SHARED_ARR) ((myComplex_t*) ((shared [] myComplex_t*) &(SHARED_ARR)[MYTHREAD]))
#define INC_MOD(NUM, INC, LIMIT) (((NUM)+(INC)) >= (LIMIT) ? ((NUM)+(INC)) - (LIMIT) : (NUM)+(INC))
#define DEC_MOD(NUM, DEC, LIMIT) (((NUM)-(DEC)) < 0 ? (LIMIT)+((NUM)-(DEC)) : (NUM)-(DEC))

static void run_1st_fft(fft3d_handle_t *handle, shared myComplex_t *dest, shared myComplex_t *source, shared myComplex_t *scratch) {
  int blk_row;
  int blk_col;
  unsigned int blk_size;
  int t,i;
  int plane;
  myComplex_t *mySource;
  myComplex_t *myScratch;
  int num_blk_rows = (handle->NY/handle->TY)/handle->blocking_factor_x;
  int num_blk_cols = handle->TY;
  
  blk_size = (handle->NX/handle->TY)*handle->blocking_factor_x;
  
  myScratch = GET_LOCAL_PORTION(scratch);
  mySource = GET_LOCAL_PORTION(source);

  MY_UPC_BARRIER();
  for(plane=0; plane<handle->num_local_planes; plane++) {
    for(blk_row=0; blk_row<num_blk_rows; blk_row++) {
      myComplex_t *currdest = myScratch+
        blk_row*blk_size*num_blk_cols+plane*num_blk_rows*num_blk_cols*blk_size;
      myComplex_t *currsrc = mySource+
        blk_row*PAD(handle->NX)*handle->blocking_factor_x+
        plane*PAD(handle->NY/handle->TY)*PAD(handle->NX);
      
      MY_UPC_POLL();
      TIMER_START(T_FFT1DROWS);
      run_1d_fft(currdest,
                 currsrc,
                 &handle->plans[0]);
      TIMER_STOP(T_FFT1DROWS);
      
      TIMER_START(T_EXCH1);
      /*make sure to reap all the old handles if we have too many of them already in flight*/
      if(handle->comm_handle_idx > COMM_HANDLE_THRESHOLD) WAIT_FOR_PUTS(handle);  
      
      if(handle->PX > 1 && handle->PY > 1 && CYCLE_EXCHANGE) {
	/*if we have a 2d processor grid then run the cycle exchange*/
	int i,j;
	for(j=1; j<=handle->PY/2; j++) {
	  int y_inc, y_dec, x_inc, x_dec;
	  y_inc = INC_MOD(handle->myloc[1], j, handle->PY);
	  y_dec = DEC_MOD(handle->myloc[1], j, handle->PY);
	  if(y_inc == y_dec) y_dec = handle->myloc[1];
	  
	  for(i=1; i<=handle->PX/2; i++) {
	    x_inc = INC_MOD(handle->myloc[0], i, handle->PX);
	    x_dec = DEC_MOD(handle->myloc[0], i, handle->PX);
	    if(x_inc == x_dec) x_dec = handle->myloc[0];
	    	  
	    int peer1 = y_inc*handle->PX + x_inc;
	    int peer2 = y_inc*handle->PX + x_dec;
	    int peer3 = y_dec*handle->PX + x_inc;
	    int peer4 = y_dec*handle->PX + x_dec;
	    
	    shared [] myComplex_t* peer_arr1;
	    shared [] myComplex_t* peer_arr2;
	    shared [] myComplex_t* peer_arr3;
	    shared [] myComplex_t* peer_arr4;
	    
	    //if(MYTHREAD==5) printf("%d sending to %d %d %d %d\n", MYTHREAD, peer1, peer2, peer3, peer4);
	    peer_arr1=
	      ((shared[] myComplex_t*) &dest[handle->mythreadplane*handle->TY + peer1])+
	      handle->mythreadrow*(handle->NY/handle->TY)*(handle->NX/handle->TY)+
	      plane*(handle->NX*(handle->NY/handle->TY));
	    
	    DO_PUT(handle, peer_arr1 + blk_row*blk_size, currdest+peer1*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	    
	    peer_arr2=
	      ((shared[] myComplex_t*) &dest[handle->mythreadplane*handle->TY + peer2])+
	      handle->mythreadrow*(handle->NY/handle->TY)*(handle->NX/handle->TY)+
	      plane*(handle->NX*(handle->NY/handle->TY));
	    
	    DO_PUT(handle, peer_arr2 + blk_row*blk_size, currdest+peer2*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	    
	    peer_arr3=
	      ((shared[] myComplex_t*) &dest[handle->mythreadplane*handle->TY + peer3])+
	      handle->mythreadrow*(handle->NY/handle->TY)*(handle->NX/handle->TY)+
	      plane*(handle->NX*(handle->NY/handle->TY));
	    
	    DO_PUT(handle, peer_arr3 + blk_row*blk_size, currdest+peer3*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	    
	    peer_arr4=
	      ((shared[] myComplex_t*) &dest[handle->mythreadplane*handle->TY + peer4])+
	      handle->mythreadrow*(handle->NY/handle->TY)*(handle->NX/handle->TY)+
	      plane*(handle->NX*(handle->NY/handle->TY));
	    
	    DO_PUT(handle, peer_arr4 + blk_row*blk_size, currdest+peer4*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	  }
	}
      } else {
	
	for(blk_col=1; blk_col<=num_blk_cols; blk_col++) {
	  int peer_thread = ((blk_col+handle->mythreadrow)%handle->TY);
	  
	  shared [] myComplex_t* peer_arr=
	    ((shared[] myComplex_t*) &dest[handle->mythreadplane*handle->TY + peer_thread])+
	    handle->mythreadrow*(handle->NY/handle->TY)*(handle->NX/handle->TY)+
	    plane*(handle->NX*(handle->NY/handle->TY));
	  
	  DO_PUT(handle, peer_arr + blk_row*blk_size, currdest+peer_thread*blk_size, 
		 sizeof(myComplex_t)*blk_size); 
	}
      }
      TIMER_STOP(T_EXCH1);
    }
  }

  TIMER_START(T_EXCH1_WAIT);
  WAIT_FOR_PUTS(handle);
  TIMER_STOP(T_EXCH1_WAIT);

  MY_UPC_BARRIER();
}

static void local_transpose_1(fft3d_handle_t *handle, shared myComplex_t *dest, 
                              shared myComplex_t *source) {

  int i, j, k, t, plane;
  unsigned int blk_size = (handle->NX/handle->TY)*(handle->blocking_factor_x);
  myComplex_t *mysrc = GET_LOCAL_PORTION(source); 
  myComplex_t *mydst = GET_LOCAL_PORTION(dest);
  
  TIMER_START(T_LOCAL_TRANSPOSE1);
  for(plane=0; plane<handle->num_local_planes; plane++) {
    for(j=0; j<handle->NX/handle->TY; j++) {
      for(k=0,i=0; i<handle->NY; i+=handle->blocking_factor_x,k++) {
        memcpy(mydst+j*PAD(handle->NY)+i, mysrc+k*blk_size+j*handle->blocking_factor_x,
               sizeof(myComplex_t)*handle->blocking_factor_x);
      }
    }
    mydst+=(PAD(handle->NX/handle->TY)*PAD(handle->NY));
    mysrc+=(handle->NX*(handle->NY/handle->TY));
  }
  TIMER_STOP(T_LOCAL_TRANSPOSE1);
}

static void run_2nd_fft(fft3d_handle_t *handle, shared myComplex_t *dest,
                        shared myComplex_t *source, shared myComplex_t *scratch) {
  int blk_row;
  int blk_col; 
  unsigned int blk_size;
  int t,plane,i;
  myComplex_t *myScratch = GET_LOCAL_PORTION(scratch);
  myComplex_t *mySource = GET_LOCAL_PORTION(source);
  
  int num_blk_rows = (handle->NX/handle->TY)/handle->blocking_factor_y;
  int num_blk_cols = handle->TZ;
  unsigned int src_plane_size = PAD(handle->NX/handle->TY)*PAD(handle->NY);
  unsigned int dst_plane_size = (handle->NY/handle->TZ)*(handle->NX/handle->TY);
  int handleidx=0;
  blk_size = (handle->NY/handle->TZ)*handle->blocking_factor_y;
  
  MY_UPC_BARRIER();
  
  for(plane=0; plane < handle->num_local_planes; plane++) {
    for(blk_row = 0; blk_row < num_blk_rows; blk_row++) {
      myComplex_t *currdest = myScratch+blk_row*blk_size*num_blk_cols+plane*src_plane_size;
      //printf("%d> %d %d\n", MYTHREAD, plane, blk_row);
      myComplex_t *currsrc = mySource+
        blk_row*PAD(handle->NY)*handle->blocking_factor_y+
        plane*src_plane_size;

      MY_UPC_POLL();
      TIMER_START(T_FFT1DCOLS);
      run_1d_fft(currdest,
                 currsrc,
                 &handle->plans[1]);
      TIMER_STOP(T_FFT1DCOLS);
      
      TIMER_START(T_EXCH2);
      /*make sure to reap all the old handles if we have too many of them already in flight*/
      if(handle->comm_handle_idx > COMM_HANDLE_THRESHOLD) WAIT_FOR_PUTS(handle);  
      
      if(handle->PZ > 1 && handle->PT > 1 && CYCLE_EXCHANGE) {
	int i,j;
	for(j=1; j<=handle->PT/2; j++) {
	  int t_inc, t_dec, z_inc, z_dec;
	  t_inc = INC_MOD(handle->myloc[3], j, handle->PT);
	  t_dec = DEC_MOD(handle->myloc[3], j, handle->PT);
	  if(t_inc==t_dec) t_dec = handle->myloc[3];
	  
	  for(i=1; i<=handle->PZ/2; i++) {
	    z_inc = INC_MOD(handle->myloc[2], i, handle->PZ);
	    z_dec = DEC_MOD(handle->myloc[2], i, handle->PZ);
	    if(z_inc==z_dec) z_dec = handle->myloc[2];
	    
	    int peer1 = t_inc*handle->PZ + z_inc;
	    int peer2 = t_inc*handle->PZ + z_dec;
	    int peer3 = t_dec*handle->PZ + z_inc;
	    int peer4 = t_dec*handle->PZ + z_dec;
	    
	    shared [] myComplex_t* peer_arr1;
	    shared [] myComplex_t* peer_arr2;
	    shared [] myComplex_t* peer_arr3;
	    shared [] myComplex_t* peer_arr4;
	    
	    peer_arr1=
	      ((shared [] myComplex_t*) &dest[peer1*handle->TY + handle->mythreadrow])+
	      handle->mythreadplane*(handle->num_local_planes)*dst_plane_size+
	      plane*dst_plane_size;
	    
	    DO_PUT(handle, peer_arr1 + blk_row*blk_size, currdest+peer1*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	    
	    peer_arr2=
	      ((shared [] myComplex_t*) &dest[peer2*handle->TY + handle->mythreadrow])+
	      handle->mythreadplane*(handle->num_local_planes)*dst_plane_size+
	      plane*dst_plane_size;
	    
	    DO_PUT(handle, peer_arr2 + blk_row*blk_size, currdest+peer2*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	    
	    peer_arr3=
	      ((shared [] myComplex_t*) &dest[peer3*handle->TY + handle->mythreadrow])+
	      handle->mythreadplane*(handle->num_local_planes)*dst_plane_size+
	      plane*dst_plane_size;
	    	    
	    DO_PUT(handle, peer_arr3 + blk_row*blk_size, currdest+peer3*blk_size, 
		   sizeof(myComplex_t)*blk_size);
	    
	    peer_arr4=
	      ((shared [] myComplex_t*) &dest[peer4*handle->TY + handle->mythreadrow])+
	      handle->mythreadplane*(handle->num_local_planes)*dst_plane_size+
	      plane*dst_plane_size;
	    
	    DO_PUT(handle, peer_arr4 + blk_row*blk_size, currdest+peer4*blk_size, 
		   sizeof(myComplex_t)*blk_size);

	  }
	}
      } else {
	
	for(blk_col = 1; blk_col <=num_blk_cols; blk_col++) {
	  int peer_thread = ((blk_col+handle->mythreadplane) % handle->TZ);
	  shared [] myComplex_t *peer_arr =
	    ((shared [] myComplex_t*) &dest[peer_thread*handle->TY + handle->mythreadrow])+
	    handle->mythreadplane*(handle->num_local_planes)*dst_plane_size+
	    plane*dst_plane_size;
	  
	  DO_PUT(handle, peer_arr + blk_row*blk_size,
		 currdest + peer_thread*blk_size, 
		 sizeof(myComplex_t)*blk_size);
	  
	}
      }
      TIMER_STOP(T_EXCH2);
    }
  }

  TIMER_START(T_EXCH2_WAIT);
  WAIT_FOR_PUTS(handle);
  TIMER_STOP(T_EXCH2_WAIT);

  MY_UPC_BARRIER();
}

static void local_transpose_2(fft3d_handle_t *handle, shared myComplex_t *dest, 
                              shared myComplex_t *source) {
  int plane;
  unsigned int dstidx;
  int blk_row, blk_col, row_idx;
  myComplex_t *mySource = GET_LOCAL_PORTION(source);
  myComplex_t *myDest = GET_LOCAL_PORTION(dest);
  unsigned int blksize = handle->blocking_factor_y;
  
  TIMER_START(T_LOCAL_TRANSPOSE2);
  for(plane=0; plane<handle->NZ; plane++) {
    
    for(dstidx=0, row_idx=0; row_idx < (handle->NY/handle->TZ); row_idx++) {
      for(blk_col = 0; blk_col < (handle->NX/handle->TY)/blksize; blk_col++, dstidx+=blksize) {
        unsigned int srcidx =
          blk_col*(handle->NY/handle->TZ)*blksize+
          row_idx*blksize;
        memcpy(myDest + dstidx, mySource+srcidx, sizeof(myComplex_t)*blksize);
      }
    }
    
    mySource += ((handle->NX/handle->TY)/blksize) * (handle->NY/handle->TZ) * blksize;
    myDest   += ((handle->NX/handle->TY)/blksize) * (handle->NY/handle->TZ) * blksize;
  }
  TIMER_STOP(T_LOCAL_TRANSPOSE2);
}

static void run_3rd_fft(fft3d_handle_t *handle, shared myComplex_t *dest, shared myComplex_t *source) {
  myComplex_t *mySource = GET_LOCAL_PORTION(source);
  myComplex_t *myDest = GET_LOCAL_PORTION(dest);
  int i;

  MY_UPC_BARRIER();
  TIMER_START(T_FFT1DPOST);
  for(i=0; i<(handle->NY/handle->TZ); i++) {
    run_1d_fft(myDest+i*PAD(handle->NX/handle->TY)*PAD(handle->NZ), 
	       mySource+i*(handle->NX/handle->TY), &handle->plans[2]);
  }
  TIMER_STOP(T_FFT1DPOST);
  MY_UPC_BARRIER();
}

void run3DFFT(fft3d_handle_t *handle, shared myComplex_t *out, shared myComplex_t *in, shared myComplex_t *scratch) {

  
  run_1st_fft(handle, out, in, scratch);
  
  local_transpose_1(handle, scratch, out);
  
  run_2nd_fft(handle, in, scratch, out);
  
  if(handle->blocking_factor_y == (handle->NX/handle->TY)) {
      run_3rd_fft(handle, out, in);
  } else {
    local_transpose_2(handle, scratch, in);
    run_3rd_fft(handle, out, scratch);
  }
}

