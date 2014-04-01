if INSTALL_OLD
bin_PROGRAMS = msprior

dist_bin_SCRIPTS = msbayes.pl

msprior_SOURCES = \
	msprior.c \
	msprior.h \
	setup.c \
	setup.h \
	../initvars.c \
	../initvars.h \
	../whiteSpaces.c \
	../whiteSpaces.h \
	../hashtab/hashtab.c \
	../hashtab/hashtab.h \
	../stringUtils.c \
	../stringUtils.h \
	../abacus/src/array_utils.c \
	../abacus/src/array_utils.h \
	../abacus/src/probability.c \
	../abacus/src/probability.h \
	../abacus/src/parsing.c \
	../abacus/src/parsing.h \
	../abacus/src/partition_combinatorics.c \
	../abacus/src/partition_combinatorics.h \
	../abacus/src/partition_combinatorics_random.c \
	../abacus/src/partition_combinatorics_random.h

msprior_LDADD = @GSL_LIBS@
msprior_CFLAGS = @GSL_CFLAGS@
endif

