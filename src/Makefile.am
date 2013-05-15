bin_PROGRAMS = \
	eureject \
	sumstatsvector \
	msDQH \
	msprior \
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

eureject_SOURCES = \
	eureject.c \
	eureject.h \
	eureject_main.c \
	array_utils.c \
	array_utils.h \
	stats_utils.c \
	stats_utils.h \
	parsing.c \
	parsing.h

sumstatsvector_SOURCES = \
	msStatsDQH.c \
	sumStatsVector.c \
	sumStatsVector.h \
	hashtab.c \
	hashtab.h \
	whiteSpaces.c \
	whiteSpaces.h \
	stringUtils.c \
	stringUtils.h

msDQH_SOURCES = \
	msQHarbpop.c \
	msQHarbpop.h \
	QHarbpop-streec.c \
	QHsubs.c \
	QHsubs.h \
	QHarbpop-getpars.c \
	rand2new.c

msprior_SOURCES = \
	msprior.c \
	msprior.h \
	setup.c \
	setup.h \
	initvars.c \
	initvars.h \
	whiteSpaces.c \
	whiteSpaces.h \
	hashtab.c \
	hashtab.h \
	stringUtils.c \
	stringUtils.h \
	array_utils.c \
	array_utils.h \
	partition_combinatorics.c \
	partition_combinatorics.h

msReject_SOURCES = rejectinC.c

# eureject_LDADD = @MATHLDADD@
# # sumstatsvector_LDADD = @LIBS@
# # msDQH_LDADD = @LIBS@
# msReject_LDADD = @MATHLDADD@
# # msprior_LDADD = @LIBS@

