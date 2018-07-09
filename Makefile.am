AUTOMAKE_OPTIONS = foreign
SUBDIRS = libsrc src 
EXTRA_DIST = data doc pfiles xspec

all-local: 
if DO_HEADAS	
	cd xspec && hmake && cd ..
else
	echo No HEADAS installed so cannot make XSPEC
endif

clean-local:
if DO_HEADAS	
	cd xspec && hmake clean && cd ..
endif

distclean-local:
if DO_HEADAS	
	cd xspec && hmake clean && cd ..
endif


