# This file is processed by `autoconf`, `autoheader`, `automake`, and `aclocal`
# to create `configure` script, `config.h.in`, `Makefile.in`s, and
# `aclocal.m4`, respectively.

# initialize Autoconf
AC_INIT([dpp-msbayes], [0.1], [joaks1@gmail.com])
AC_PREREQ([2.59])

# safety check for source dir
AC_CONFIG_SRCDIR([src/msprior.c])

# designate dir for auxiliary scripts
AC_CONFIG_AUX_DIR([build-aux])

# set default install prefix
AC_PREFIX_DEFAULT([/usr/local])

# initialize Automake
AM_INIT_AUTOMAKE([1.9 foreign -Wall dist-tarZ tar-ustar filename-length-max=299])

# check for C compiler
AC_PROG_CC
# support per-target flags
AM_PROG_CC_C_O

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

AC_ARG_ENABLE([static],
    [AS_HELP_STRING(
        [--enable-static],
        [build statically linked binaries])],
    [STATIC_ARG="${enable_static}"],
    [STATIC_ARG=no])

AC_ARG_ENABLE([old-model],
    [AS_HELP_STRING(
        [--enable-old-model],
        [build and install old model along side new implementation])],
    [OLD_MODEL_ARG="${enable_old_model}"],
    [OLD_MODEL_ARG=no])

if ! test "${STATIC_ARG}" = "no"
then
    LDFLAGS="${LDFLAGS} -static"
fi

AM_CONDITIONAL([INSTALL_OLD], [test "$OLD_MODEL_ARG" = "yes"])

###############################################################################
# check for GSL
###############################################################################

AC_ARG_WITH([gsl],
    [AS_HELP_STRING(
        [--with-gsl=yes/no/DIR],
        [if yes, configure will look for a system-level GSL installation (this is the default); if no, configure will only build and install gsl-independent ABACUS tools; if a directory is provided, configure will use the gsl library installed at that location. For example, if the non-system-level GSL library files you wish to use are in /home/joebob/lib, then use --with-gsl=/home/joebob.])],
    [GSL_ARG="${with_gsl}"],
    [GSL_ARG=yes])

GSL_FOUND="no"
GSL_LIBS=
GSL_CFLAGS=
if test "${GSL_ARG}" = "yes"
then
    GSL_CONFIG="gsl-config"
    AC_MSG_NOTICE([looking for gsl using $GSL_CONFIG...])
    if test -z "`${GSL_CONFIG} --version || echo ''`"
    then
        AC_MSG_NOTICE([could not find GSL library with $GSL_CONFIG])
    else
        GSL_LIBS="`${GSL_CONFIG} --libs`"
        GSL_CFLAGS="`${GSL_CONFIG} --cflags`"
        GSL_FOUND="yes"
        AC_MSG_NOTICE([GSL library found with $GSL_CONFIG])
    fi

    if test "$GSL_FOUND" = "no"
    then
        GSL_CONFIG="/usr/bin/gsl-config"
        AC_MSG_NOTICE([looking for gsl using $GSL_CONFIG...])
        if test -z "`${GSL_CONFIG} --version || echo ''`"
        then
            AC_MSG_NOTICE([could not find GSL library with $GSL_CONFIG])
        else
            GSL_LIBS="`${GSL_CONFIG} --libs`"
            GSL_CFLAGS="`${GSL_CONFIG} --cflags`"
            GSL_FOUND="yes"
            AC_MSG_NOTICE([GSL library found with $GSL_CONFIG])
        fi
    fi

    if test "$GSL_FOUND" = "no"
    then
        GSL_CONFIG="/usr/local/bin/gsl-config"
        AC_MSG_NOTICE([looking for gsl using $GSL_CONFIG...])
        if test -z "`${GSL_CONFIG} --version || echo ''`"
        then
            AC_MSG_NOTICE([could not find GSL library with $GSL_CONFIG])
        else
            GSL_LIBS="`${GSL_CONFIG} --libs`"
            GSL_CFLAGS="`${GSL_CONFIG} --cflags`"
            GSL_FOUND="yes"
            AC_MSG_NOTICE([GSL library found with $GSL_CONFIG])
        fi
    fi

    if test "$GSL_FOUND" = "no"
    then
        AC_MSG_NOTICE([looking for GSL using pkg-config...])
        if test -z "${PKG_CONFIG}"
        then
            AC_PATH_PROG(PKG_CONFIG, pkg-config)
        fi
        if test -n "${PKG_CONFIG}"
        then
            GSL_CONFIG="${PKG_CONFIG} libgsl"
            if test -z "`${GSL_CONFIG} --modversion || echo ''`"
            then
                AC_MSG_NOTICE([could not find GSL library with pkg-config])
            else
                GSL_LIBS="`${GSL_CONFIG} --libs`"
                GSL_CFLAGS="`${GSL_CONFIG} --cflags`"
                GSL_FOUND="yes"
                AC_MSG_NOTICE([GSL library found with pkg-config])
            fi
        else
            AC_MSG_NOTICE([could not find pkg-config to search for GSL])
        fi
    fi

    if test "$GSL_FOUND" = "no"
    then
        AC_MSG_ERROR([could not find GSL library. Either install GSL or designate its location with --with-gsl=path/to/gsl/base/dir. If you do not want to deal with GSL you can also build and install only GSL-independent tools using the --without-gsl configure flag.])
    fi

elif test "${GSL_ARG}" = "no"
then
    AC_MSG_NOTICE([--without-gsl specified; configure will ONLY build and install GSL-independent ABACUS tools.])
else
    if ! test -d "${GSL_ARG}/lib"
    then
        AC_MSG_ERROR([Expecting a directory named $GSL_ARG/lib; please verify --with-gsl argument is correct.])
    else
        GSL_CONFIG="${GSL_ARG}/bin/gsl-config"
        AC_MSG_NOTICE([looking for $GSL_CONFIG...])
        if test -z "`${GSL_CONFIG} --version || echo ''`"
        then
            AC_MSG_NOTICE([could not find gsl-config at --with-gsl specified location. Trying basic GSL configuration based on location])
            GSL_LIBS="-L${GSL_ARG}/lib -lgsl -lgslcblas -lm"
            GSL_CFLAGS="-I${GSL_ARG}/include"
            GSL_FOUND="yes"
        else
            GSL_LIBS="`${GSL_CONFIG} --libs`"
            GSL_CFLAGS="`${GSL_CONFIG} --cflags`"
            GSL_FOUND="yes"
            AC_MSG_NOTICE([GSL library found with $GSL_CONFIG])
        fi
    fi
fi

AM_CONDITIONAL([USE_GSL], [test "$GSL_FOUND" = "yes"])
AC_SUBST([GSL_LIBS])
AC_SUBST([GSL_CFLAGS])
          
# define flags
MATHLDADD="-lm"
AC_SUBST([MATHLDADD])

AC_SUBST([LDFLAGS])

# declare `config.h` as output header
AC_CONFIG_HEADERS([config.h])

# declare output `Makefile`s
AC_CONFIG_FILES([Makefile src/Makefile src/old-model/Makefile])

# call 'abacus's configure script recursively
AC_CONFIG_SUBDIRS([src/abacus])

# output declared files
AC_OUTPUT

