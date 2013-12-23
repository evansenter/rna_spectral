CCFLAGS           = -c -std=c99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q
LDFLAGS           = -lm -lgomp -llapack -L/usr/local/include -llapacke -lgslcblas -lgsl -o
BINDIR            = ~/bin # Change this to the BINDIR
LIBDIR            = ~/lib # Change this to the LIBDIR
CC                = gcc-4.9
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

all: RNAspectral
	
RNAspectral: spectral_grid.o spectral_params.o spectral_functions.o
	$(CC) spectral_grid.o spectral_params.o spectral_functions.o -lRNA $(LDFLAGS) RNAspectral
	ar cr ../shared/libspectral.a spectral_functions.o spectral_params.o
	
spectral_grid.o: spectral_grid.c spectral_grid.h spectral_functions.h spectral_params.h ../shared/constants.h
	$(CC) $(CCFLAGS) spectral_grid.c
	
spectral_functions.o: spectral_functions.c spectral_functions.h spectral_params.h spectral_grid.h
	$(CC) $(CCFLAGS) spectral_functions.c
	
spectral_params.o: spectral_params.c spectral_params.h energy_const.h
	$(CC) $(CCFLAGS) spectral_params.c

clean:
	rm -f *.o RNAspectral
	
install: RNAspectral
	cp RNAspectral $(BINDIR)
	cp ../shared/libspectral.a $(LIBDIR)
		