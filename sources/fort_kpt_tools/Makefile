include ../make.inc

all : fort_kpt_tools.so #gaunt.so gutils.so fort_kpt_tools.so

gaunt.so: gaunt.f90
	$(CMP) -c $? -m gaunt 
	mv gaunt.*.so gaunt.so

gutils.so: gutils.cxx
	python2 setup.py build_ext --inplace

fort_kpt_tools.so: fort_kpt_tools.f
	$(CMP) -c $? -m fort_kpt_tools
	mv fort_kpt_tools.*.so fort_kpt_tools.so

clean:
	rm -rf build fort_kpt_tools.so 
	rm -rf build gutils.so gaunt.so
	rm -rf *.dSYM 

