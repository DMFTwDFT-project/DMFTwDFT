include ../make.inc

all : gaunt.so gutils.so fort_kpt_tools.so

gaunt.so: gaunt.f90
	$(CMP) -c $? -m gaunt 

gutils.so: gutils.cxx
	python2 setup.py build_ext --inplace

fort_kpt_tools.so: fort_kpt_tools.f
	$(CMP) -c $? -m fort_kpt_tools

clean:
	rm -rf build fort_kpt_tools.so 
	rm -rf build gutils.so gaunt.so
	rm -rf *.dSYM 

