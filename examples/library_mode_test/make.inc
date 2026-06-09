
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpiifort
# Preserve f2py shared-library linker flags under conda/numpy.distutils.
CMP = env CC=gcc CXX=g++ NPY_DISTUTILS_APPEND_FLAGS=1 f2py --fcompiler=gnu95 # fortran to python converter


FCOPTS=-O2 #-check 
LDOPTS=-O2 

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = $(MKLROOT)/lib/intel64
LIBS   =  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

#LIBDIR = /opt/intel/mkl/lib/intel64
#LIBS   =  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

#=======================
# ATLAS Blas and LAPACK
#=======================
#LIBDIR = /usr/local/lib
#LIBS = -L$(LIBDIR)  -llapack -lf77blas -lcblas -latlas
