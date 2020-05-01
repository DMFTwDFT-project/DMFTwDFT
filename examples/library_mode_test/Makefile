include make.inc


all : test.x
test.x : test.F90
	$(MPIF90) $(FCOPTS) -o $@ test.F90 $(LIBS) -L./ -ldmft

clean :
	rm -f test.x 


