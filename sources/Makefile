ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../Makefile.2
FORTKPTTOOLSMAKEFILE=Makefile
DMFTKSUMMAKEFILE=Makefile

all: lib DMFTwDFT dos copy

default: all

DMFTwDFT: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) dmft)
	(cd $(ROOTDIR)/fort_kpt_tools && $(MAKE) -f $(FORTKPTTOOLSMAKEFILE) fort_kpt_tools.so )
	(cd $(ROOTDIR)/dmft_ksum && $(MAKE) -f $(DMFTKSUMMAKEFILE) all )

dos: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) dos)
lib: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) libs)

copy:
	cp  *.x ../bin/
	cp  dmft_ksum/dmft_ksum_band ../bin/
	cp  dmft_ksum/dmft_ksum_partial_band ../bin/
	cp  fort_kpt_tools/fort_kpt_tools.so ../bin/

clean:
	rm -f *.x
	rm -f libdmft.a
	rm -f $(ROOTDIR)/fort_kpt_tools/fort_kpt_tools.so
	rm -f $(ROOTDIR)/dmft_ksum/dmft_ksum_band
	rm -f $(ROOTDIR)/dmft_ksum/dmft_ksum_partial_band
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f src/*~
	@( cd $(ROOTDIR) && if [ -d src/obj ] ; \
		then cd src/obj && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf obj ; \
	fi )

objdir:
	@( cd $(ROOTDIR) && if [ ! -d src/obj ] ; \
                then mkdir src/obj ; \
        fi ) ;
