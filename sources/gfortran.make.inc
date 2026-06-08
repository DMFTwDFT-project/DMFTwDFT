#=============================
# For Linux with gnu compilers 
#=============================
F90 = gfortran
COMMS=mpi
MPIF90=mpif90
# Preserve f2py shared-library linker flags under conda/numpy.distutils.
CMP = env NPY_DISTUTILS_APPEND_FLAGS=1 f2py --fcompiler=gnu95 # fortran to python converter

FCOPTS=-O2 -ffree-line-length-none #-check
LDOPTS=-O2

#=======================
# Blas and LAPACK
#=======================

LIBDIR = /usr/local/lib/ 
LIBS = -L$(LIBDIR)  -llapack -lblas
