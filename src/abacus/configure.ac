# This file is processed by `autoconf`, `autoheader`, `automake`, and `aclocal`
# to create `configure` script, `config.h.in`, `Makefile.in`s, and
# `aclocal.m4`, respectively.

# initialize Autoconf
AC_INIT([abacus], [0.1], [joaks1@ku.edu])
AC_PREREQ([2.59])

# safety check for source dir
AC_CONFIG_SRCDIR([src/eureject.c])

# designate dir for auxiliary scripts
AC_CONFIG_AUX_DIR([build-aux])

# set default install prefix
AC_PREFIX_DEFAULT([/usr/local])

# initialize Automake
AM_INIT_AUTOMAKE([1.9 foreign -Wall dist-tarZ tar-ustar filename-length-max=299])

# check for C compiler
AC_PROG_CC

# check for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h stdlib.h math.h assert.h unistd.h ctype.h string.h]) 

# check for library functions
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([pow strtok strncpy calloc sqrt]) 

# check for "Check"
# AM_PATH_CHECK([0.9.0],[],[AC_MSG_WARN([check is missing or old; you will not be able to run the test suite])])
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [], [AC_MSG_WARN([check is missing or old; you will not be able to run the test suite])])

###############################################################################
# check for GSL
# (from http://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html)
# These tests create macros:
#   HAVE_LIBGSL, HAVE_LIBGSLCBGLAS, HAVE_LIBM
# They will also add '-lgls -lgslcblas -lm' to `LIBS`
###############################################################################
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
          
# define flags
MATHLDADD="-lm"
AC_SUBST([MATHLDADD])

# declare `config.h` as output header
AC_CONFIG_HEADERS([config.h])

# declare output `Makefile`s
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])

# output declared files
AC_OUTPUT

