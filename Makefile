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
GCC_GTEQ_4.9.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40900)
LIB              := ../../lib
H                := ../../h

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	CCFLAGS += -fdiagnostics-color=always
endif
	
RNAspectral.out: rna_spectral.o spectral_params.o spectral_functions.o spectral_initializers.o
	$(CC) rna_spectral.o spectral_params.o spectral_functions.o spectral_initializers.o -lRNA $(LDFLAGS) RNAspectral.out
	ar cr $(LIB)/libspectral.a spectral_functions.o spectral_params.o spectral_initializers.o
	
rna_spectral.o: rna_spectral.c $(H)/vienna_functions.h $(H)/spectral_functions.h $(H)/spectral_params.h $(H)/constants.h $(H)/spectral_initializers.h
	$(CC) $(CCFLAGS) rna_spectral.c
	
spectral_functions.o: spectral_functions.c $(H)/spectral_functions.h $(H)/spectral_params.h $(H)/spectral_data_structures.h $(H)/spectral_initializers.h
	$(CC) $(CCFLAGS) spectral_functions.c
	
spectral_params.o: spectral_params.c $(H)/spectral_params.h $(H)/spectral_data_structures.h
	$(CC) $(CCFLAGS) spectral_params.c
  
spectral_initializers.o: spectral_initializers.c $(H)/spectral_initializers.h $(H)/spectral_functions.h
	$(CC) $(CCFLAGS) spectral_initializers.c
  
$(H)/spectral_data_structures.h: $(H)/vienna_data_structures.h $(H)/energy_const.h

clean:
	rm -f *.o $(LIB)/libspectral.a RNAspectral.out
	
install: RNAspectral.out
	cp RNAspectral.out $(BINDIR)/RNAspectral
	cp $(LIB)/libspectral.a $(LIBDIR)
		
