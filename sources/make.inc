#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpif90
CMP = f2py --fcompiler=intelem # fortran to python converter

FCOPTS=-O2 -fallow-argument-mismatch -ffree-line-length-none -L /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem #-check 
LDOPTS=-O2 

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = $(MKLROOT)/lib/intel64
LIBS   =  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

