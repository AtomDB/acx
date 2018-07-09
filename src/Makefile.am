AM_CPPFLAGS = -g -I. -I../libsrc/ape/include -I../libsrc/cfitsio \
	-DDIRECTORY='"@abs_top_builddir@"'

bin_PROGRAMS = dacx

dacx_SOURCES = dacx.c get_input_data.c messages.c expand_env.c get_abund.c \
	calc_cx.c messages.h dacx.h acx.h

dacx_LDADD = ../libsrc/cfitsio/libcfitsio.a ../libsrc/ape/src/libape.a

install-exec-hook:
#	cp -p dacx ../bin/
	mkdir -p @prefix@/syspfiles
	cp -p ../pfiles/dacx.par @prefix@/syspfiles/
