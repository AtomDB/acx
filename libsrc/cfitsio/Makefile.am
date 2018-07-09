#
#  Revised and bastardized version by Randall Smith.  Sorry.
#
#
# Makefile for cfitsio library:
#       libcfits.a
#
# Oct-96 : original version by 
#
#       JDD/WDP
#       NASA GSFC
#       Oct 1996
#
# 25-Jan-01 : removed conditional drvrsmem.c compilation because this
#             is now handled within the source file itself.
# 09-Mar-98 : modified to conditionally compile drvrsmem.c. Also
# changes to target all (deleted clean), added DEFS, LIBS, added
# DEFS to .c.o, added SOURCES_SHMEM and MY_SHMEM, expanded getcol*
# and putcol* in SOURCES, modified OBJECTS, mv changed to /bin/mv
# (to bypass aliasing), cp changed to /bin/cp, add smem and
# testprog targets. See also changes and comments in configure.in
#

SHELL =		/bin/sh
RANLIB =	@RANLIB@
CC =		@CC@
SSE_FLAGS =	@SSE_FLAGS@
FC =		@FC@
DEFS =		@DEFS@
LIBS =		@LIBS@
FLEX =		flex
BISON =		bison

SHLIB_LD =	@SHLIB_LD@
SHLIB_SUFFIX =	@SHLIB_SUFFIX@

lib_LIBRARIES = libcfitsio.a

libcfitsio_a_SOURCES = 	buffers.c cfileio.c checksum.c drvrfile.c drvrmem.c \
		drvrnet.c drvrsmem.c drvrgsiftp.c editcol.c edithdu.c eval_l.c \
		eval_y.c eval_f.c fitscore.c getcol.c getcolb.c getcold.c getcole.c \
		getcoli.c getcolj.c getcolk.c getcoll.c getcols.c getcolsb.c \
		getcoluk.c getcolui.c getcoluj.c getkey.c group.c grparser.c \
		histo.c iraffits.c \
		modkey.c putcol.c putcolb.c putcold.c putcole.c putcoli.c \
		putcolj.c putcolk.c putcoluk.c putcoll.c putcols.c putcolsb.c \
		putcolu.c putcolui.c putcoluj.c putkey.c region.c scalnull.c \
		swapproc.c wcssub.c wcsutil.c imcompress.c quantize.c ricecomp.c \
		pliocomp.c fits_hcompress.c fits_hdecompress.c zuncompress.c zcompress.c \
		adler32.c crc32.c inffast.c inftrees.c trees.c zutil.c \
		deflate.c infback.c inflate.c uncompr.c \
		cfortran.h    drvrsmem.h   fitsio.h   grparser.h  inftrees.h  zconf.h \
		crc32.h       eval_defs.h  fitsio2.h  inffast.h   longnam.h   zlib.h \
		deflate.h     eval_tab.h   fpack.h    inffixed.h  region.h    zutil.h \
		drvrgsiftp.h  f77_wrap.h   group.h    inflate.h   trees.h 

EXTRA_DIST =	f77_wrap1.c f77_wrap2.c f77_wrap3.c f77_wrap4.c f77.inc

bin_PROGRAMS = smem testprog fpack funpack fitscopy speed imcopy listhead cookbook

smem_SOURCES = smem.c
smem_LDADD = libcfitsio.a

testprog_SOURCES = testprog.c
testprog_LDADD = libcfitsio.a

fpack_SOURCES = fpack.c fpackutil.c
fpack_LDADD = libcfitsio.a

funpack_SOURCES = funpack.c fpackutil.c
funpack_LDADD = libcfitsio.a

fitscopy_SOURCES = fitscopy.c
fitscopy_LDADD = libcfitsio.a

speed_SOURCES = speed.c
speed_LDADD = libcfitsio.a

imcopy_SOURCES = imcopy.c
imcopy_LDADD = libcfitsio.a

listhead_SOURCES = listhead.c
listhead_LDADD = libcfitsio.a

cookbook_SOURCES = cookbook.c
cookbook_LDADD = libcfitsio.a

eval:		# Rebuild eval_* files from flex/bison source
		$(FLEX) -t eval.l > eval_l.c1
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' eval_l.c1 > eval_l.c
		/bin/rm -f eval_l.c1
		$(BISON) -d -v -y eval.y
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' y.tab.c > eval_y.c
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' y.tab.h > eval_tab.h
		/bin/rm -f y.tab.c y.tab.h

