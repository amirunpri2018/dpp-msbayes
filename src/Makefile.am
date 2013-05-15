bin_PROGRAMS = eureject

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

eureject_LDADD = $(MATHLDADD)

