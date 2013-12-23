# Makefile for RNAspectral

CCFLAGS           = -c -std=c99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I ../../h
LDFLAGS           = -lm -lgomp -llapack -L/usr/local/include -llapacke -lgslcblas -lgsl -o
BINDIR            = ~/bin # Change this to the BINDIR
LIBDIR            = ~/lib # Change this to the LIBDIR
CC                = gcc
GCC_VERSION      := $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  := $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40600)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif
	
RNAspectral.out: spectral_grid.o spectral_params.o spectral_functions.o
	$(CC) spectral_grid.o spectral_params.o spectral_functions.o -lRNA $(LDFLAGS) RNAspectral.out
	ar cr ../../lib/libspectral.a spectral_functions.o spectral_params.o
	
spectral_grid.o: spectral_grid.c ../../h/spectral_grid.h ../../h/spectral_functions.h ../../h/spectral_params.h ../../h/constants.h
	$(CC) $(CCFLAGS) spectral_grid.c
	
spectral_functions.o: spectral_functions.c ../../h/spectral_functions.h ../../h/spectral_params.h ../../h/spectral_grid.h
	$(CC) $(CCFLAGS) spectral_functions.c
	
spectral_params.o: spectral_params.c ../../h/spectral_params.h ../../h/energy_const.h
	$(CC) $(CCFLAGS) spectral_params.c

clean:
	rm -f *.o ../../lib/libspectral.a RNAspectral.out
	
install: RNAspectral.out
	cp RNAspectral.out $(BINDIR)/RNAspectral
	cp ../../lib/libspectral.a $(LIBDIR)
		
