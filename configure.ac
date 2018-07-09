#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([acx], [1.0.1], [rsmith@cfa.harvard.edu])
AC_CONFIG_SRCDIR([src/dacx.c])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow strchr])

# Checks for libraries, builds them if need be.
AC_CHECK_LIB([ape], [check_args],
		    [],
		    [AC_CONFIG_SUBDIRS([libsrc/ape])])
#AC_CHECK_LIB([cfitsio], [ffgpxv],
#			[],
#			[AC_CONFIG_SUBDIRS([libsrc/cfitsio])])

AC_CONFIG_SUBDIRS([libsrc/cfitsio])

AC_CONFIG_FILES([xspec/acx.h])
AC_CONFIG_FILES([libsrc/Makefile
                 src/Makefile
		 Makefile])
AS_IF([test "x$XANBIN" != x], 
	    [AC_CONFIG_COMMANDS([xspec/lpack_acx.cxx], 
	    [cd xspec && initpackage acx model.dat . && cd ..])],
	    [echo No HEADAS installed so no XSPEC modules compiled.])
AM_CONDITIONAL([DO_HEADAS], [test "x$XANBIN" != x])

dnl AS_IF([test "x$HEADAS" != x], 
dnl	    [AC_CONFIG_COMMANDS([xspec/libacx.so],
dnl	    [cd xspec && hmake && cd ..])], 
dnl	    [echo No HEADAS installed so no XSPEC modules compiled.])

AC_OUTPUT
