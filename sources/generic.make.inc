
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = 
COMMS=
MPIF90=
# Preserve f2py shared-library linker flags under conda/numpy.distutils.
CMP = env NPY_DISTUTILS_APPEND_FLAGS=1 f2py --fcompiler=gnu95 # fortran to python converter

FCOPTS=-O2 -ffree-line-length-none #-check
LDOPTS=-O2

#=======================
# Blas and LAPACK
#=======================

LIBDIR = 
LIBS = 
