.SUFFIXES: .c .upc .o

PLATFORM ?= fftw3
PAR_BUILD ?= 0
UPCC ?= cc  
CLASS ?= DD16

include Make.$(PLATFORM)

COMM_OBJS = verify.o init_exp.o c_randdp.o myfft_$(FFT_BACKEND).o 
CFLAGS = -I. $(FFT_INCFLAGS) 
UPCFLAGS = -h upc -I. -DCLASS=$(CLASS) 
LIBS =  $(FFT_LIBS) -lm
ALL_CLASSES ?= DDE DDQ DDH DD DD2 DD4 DD8

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

.upc.o: fft3d.uph
	$(UPCC) $(UPCFLAGS) -c $< -o $@

all:
	$(MAKE) clean
	$(MAKE) upc-bench

all-bench:
	for c in $(ALL_CLASSES); do \
		$(MAKE) upc-bench CLASS=$$c;\
	done

upc-bench: $(COMM_OBJS) fft3d.upc timers.upc fft3d.uph 
	$(UPCC) $(UPCFLAGS) -c fft3d.upc -o fft3d.o
	$(UPCC) $(UPCFLAGS) -c timers.upc -o timers.o
	$(UPCC) $(UPCFLAGS) $(COMM_OBJS) fft3d.o timers.o ft-2d.upc -o ft-2d-upc.$(FFT_BACKEND).$(CLASS) $(LIBS)

clean:
	rm -rf *~ *.o *.B *.t *.i *_temps *.N *.trans.c  
