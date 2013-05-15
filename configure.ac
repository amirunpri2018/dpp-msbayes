# This file is processed by `autoconf`, `autoheader`, `automake`, and `aclocal`
# to create `configure` script, `config.h.in`, `Makefile.in`s, and
# `aclocal.m4`, respectively.

# initialize Autoconf
AC_INIT([msbayes], [0.1], [joaks1@ku.edu])
AC_PREREQ([2.59])

# safety check for source dir
AC_CONFIG_SRCDIR([src/eureject.c])

# designate dir for auxiliary scripts
AC_CONFIG_AUX_DIR([build-aux])

# initialize Automake
AM_INIT_AUTOMAKE([1.9 foreign -Wall dist-tarZ tar-ustar filename-length-max=299])

# check for C compiler
AC_PROG_CC

# declare `config.h` as output header
AC_CONFIG_HEADERS([config.h])

# declare output `Makefile`s
AC_CONFIG_FILES([Makefile src/Makefile])

# output declared files
AC_OUTPUT


