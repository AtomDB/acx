AM_CPPFLAGS = -I. -I../include -I../include/ape

bin_PROGRAMS = pcheck pget plist pquery2 pset punlearn 

lib_LIBRARIES = libape.a

libape_a_SOURCES = ape_binary.c ape_error.c ape_io.c ape_list.c ape_msg.c \
		  ape_par.c ape_par_group.c ape_test.c ape_trad.c ape_util.c \
		pil.c pil_error.c

pcheck_SOURCES = pcheck.c
pcheck_LDADD = libape.a

pget_SOURCES = pget.c
pget_LDADD = libape.a

plist_SOURCES = plist.c
plist_LDADD = libape.a

pquery2_SOURCES = pquery2.c
pquery2_LDADD = libape.a

pset_SOURCES = pset.c
pset_LDADD = libape.a

punlearn_SOURCES = punlearn.c
punlearn_LDADD = libape.a

install_exec_hook:
	cp -p pcheck ../bin
	cp -p pget ../bin
	cp -p plist ../bin
	cp -p pquery2 ../bin
	cp -p pset ../bin
	cp -p punlearn ../bin
	cp -p libape.a ../lib
