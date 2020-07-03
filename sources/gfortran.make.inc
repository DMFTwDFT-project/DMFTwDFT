
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = gfortran
COMMS=mpi
MPIF90=mpif90
CMP = f2py --fcompiler=intelem # fortran to python converter

FCOPTS=-O2 -ffree-line-length-none #-check
LDOPTS=-O2

#=======================
# Blas and LAPACK
#=======================

LIBDIR = /home/uthpala/lib
LIBS = -L$(LIBDIR)  -llapack -lblas
