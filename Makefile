ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../Makefile.2


default: DMFTwDFT 


DMFTwDFT: objdir 
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) dmft) 

dos: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) dos)
lib: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) libs)

clean:
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

