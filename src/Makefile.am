EXTRA_DIST = hashtab/main.c hashtab/Makefile

bin_PROGRAMS = \
	sumstatsvector \
	msDQH \
	dpp-msprior \
	msReject

dist_bin_SCRIPTS = \
	acceptRej.pl \
	acceptRej_flex.pl \
	convertIM.pl \
	dataSummary.pl \
	msbayes.pl \
	msCombModels.pl \
	obsSumStats.pl \
	regress_cli.r

rlibdir = $(libdir)/msbayes
dist_rlib_SCRIPTS = \
	acceptRej_flex.r \
	acceptRej.r \
	calmod.r \
	loc2plot.r \
	make_pd2005.r

sumstatsvector_SOURCES = \
	msStatsDQH.c \
	sumStatsVector.c \
	sumStatsVector.h \
	hashtab/hashtab.c \
	hashtab/hashtab.h \
	whiteSpaces.c \
	whiteSpaces.h \
	stringUtils.c \
	stringUtils.h

sumstatsvector_LDADD = @GSL_LIBS@
sumstatsvector_CFLAGS = @GSL_CFLAGS@

msDQH_SOURCES = \
	msQHarbpop.c \
	msQHarbpop.h \
	QHarbpop-streec.c \
	QHsubs.c \
	QHsubs.h \
	QHarbpop-getpars.c \
	rand2new.c

msDQH_LDADD = @GSL_LIBS@
msDQH_CFLAGS = @GSL_CFLAGS@

dpp_msprior_SOURCES = \
	msprior.c \
	msprior.h \
	setup.c \
	setup.h \
	initvars.c \
	initvars.h \
	whiteSpaces.c \
	whiteSpaces.h \
	hashtab/hashtab.c \
	hashtab/hashtab.h \
	stringUtils.c \
	stringUtils.h \
	abacus/src/array_utils.c \
	abacus/src/array_utils.h \
	abacus/src/probability.c \
	abacus/src/probability.h \
	abacus/src/parsing.c \
	abacus/src/parsing.h \
	abacus/src/partition_combinatorics.c \
	abacus/src/partition_combinatorics.h \
	abacus/src/partition_combinatorics_random.c \
	abacus/src/partition_combinatorics_random.h

dpp_msprior_LDADD = @GSL_LIBS@
dpp_msprior_CFLAGS = @GSL_CFLAGS@

msReject_SOURCES = rejectinC.c

msReject_LDADD = @MATHLDADD@

