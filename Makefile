EXECS   := rockstar_linker kali_linker 

OBJS   :=	./io.o	\
			./io_hdf5.o \
			./particles.o 

EXECS_OBJS := rockstar_main.o \
			  kali_main.o

INCL   :=	./io.h \
			./io_hdf5.h \
			./particles.h \
			./Makefile
 
USE_MPI=TRUE
#USE_MPI=FALSE

#USE_OPENMP=TRUE
USE_OPENMP=FALSE


GSL_DIR := $(shell gsl-config --prefix)
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
GSL_LIBDIR := $(GSL_DIR)/lib

HDF5INCL := -I/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/include
HDF5LIB := -L/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/lib

OPTIMIZE = -g -O0 -Wall # optimization and warning flags

LIBS   += -g -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CCFLAGS := 

ifeq ($(USE_MPI),TRUE)
    OPTS += -DMPI #  This creates an MPI version that can be used to process files in parallel
    CC := mpicc  # sets the C-compiler
	LIBS := -lmpi
else
    CC = cc  # sets the C-compiler
endif

ifeq ($(USE_OPENMP),TRUE)
    CCFLAGS += -fopenmp 
    OPTS += -DOPENMP 
endif

OPTS += -DUSE_HDF5 -DADJUST_PROCESSORS 
LIBS += $(HDF5LIB) -lhdf5
CCFLAGS += $(HDF5INCL)  


CFLAGS =   $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL) $(CCFLAGS) $(OPTS) 

all: $(EXECS)
	@if [ "$(USE_MPI)" = "TRUE" ]; then echo "RUNNING WITH MPI"; else echo "MPI DISABLED"; fi
	@if [ "$(USE_OPENMP)" = "TRUE" ]; then echo "RUNNING WITH OPENMP"; else echo "OPENMP DISABLED"; fi

rockstar_linker: $(OBJS) rockstar_main.o
	$(CC) $(CCFLAGS) $^ $(LIBS) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) -o  $@ 

kali_linker: $(OBJS) kali_main.o
	$(CC) $(CCFLAGS) $^ $(LIBS) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) -o  $@ 
clean:
	rm -f $(OBJS) $(EXECS_OBJS) $(EXECS)

.PHONY: all clean clena celan celna

celan celna clena claen:clean
	
