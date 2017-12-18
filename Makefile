EXECS   := particle_linker 

OBJS   := 	./main.o \
			./io.o	\
			./io_hdf5.o \
			./particles.o 

INCL   :=	./io.h \
			./io_hdf5.h \
			./particles.h \
			./Makefile
 
#USE_MPI=TRUE
USE_MPI=FALSE

ifeq ($(USE_MPI),TRUE)
    OPTS += -DMPI #  This creates an MPI version that can be used to process files in parallel
    CC := mpicc  # sets the C-compiler
	LIBS := -lmpi
else
    CC = cc  # sets the C-compiler
endif


OPTS += -DUSE_HDF5 
LIBS += $(HDF5LIB) -lhdf5
CFLAGS += $(HDF5INCL)  

GSL_DIR := $(shell gsl-config --prefix)
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
GSL_LIBDIR := $(GSL_DIR)/lib

HDF5INCL := -I/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/include
HDF5LIB := -L/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/lib

OPTIMIZE = -g -O0 -Wall -Werror # optimization and warning flags

OPTS += 

LIBS   += -g -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS =   $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL) $(OPTS) 

all: $(EXECS)
	@if [ "$(USE_MPI)" = "TRUE" ]; then echo "RUNNING WITH MPI"; else echo "MPI DISABLED"; fi

particle_linker: $(OBJS)
	$(CC) $(CCFLAGS) $^ $(LIBS) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) -o  $@ 

clean:
	rm -f $(OBJS) $(EXECS)

.PHONY: all clean clena celan celna

celan celna clena claen:clean
	
