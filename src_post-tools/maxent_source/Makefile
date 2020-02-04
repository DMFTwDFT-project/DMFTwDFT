include ../../make.inc
# CMP = f2py #--fcompiler=intelem
#CMPLIBS = --opt='-fast' --link-lapack_opt
#CMPLIBS2 = --f90flags='-openmp ' --opt='-fast' --link-lapack_opt

all : maxent_routines.so maxent_routines_OMP.so

maxent_routines.so : maxent_routines.f90
	$(CMP) -c maxent_routines.f90 -m maxent_routines $(CMPLIBS)

maxent_routines_OMP.so : maxent_routines_OMP.f90
	$(CC) -c cputim.c -fPIC
	$(CMP) -c maxent_routines_OMP.f90 -m maxent_routines_OMP cputim.o $(CMPLIBS2)

clean :
	rm -f maxent_routines.so maxent_routines_OMP.so *.pyc cputim.o
	rm -rf *.dSYM


