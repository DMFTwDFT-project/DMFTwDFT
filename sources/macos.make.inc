#===============================================
# For macOS x86_64 conda with gnu compilers
#===============================================
F90 = gfortran
COMMS=mpi
MPIF90=mpif90
# Preserve f2py shared-library linker flags under conda/numpy.distutils.
CMP = env SDKROOT=$(shell xcrun --show-sdk-path) CC=/usr/bin/clang CXX=/usr/bin/clang++ F77=gfortran F90=gfortran FC=gfortran FFLAGS="-arch x86_64 -fPIC -O3 -funroll-loops" CFLAGS="-arch x86_64 -isysroot $(shell xcrun --show-sdk-path)" LDFLAGS="-undefined dynamic_lookup -bundle -arch x86_64 -Wl,-syslibroot,$(shell xcrun --show-sdk-path)" f2py --fcompiler=gnu95 # fortran to python converter

FCOPTS=-O2 -ffree-line-length-none #-check
LDOPTS=-O2 -Wl,-syslibroot,$(shell xcrun --show-sdk-path)

#=======================
# Blas and LAPACK
#=======================

LIBDIR = /usr/local/lib/
LIBS = -L$(LIBDIR)  -llapack -lblas
