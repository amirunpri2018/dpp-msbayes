dist_noinst_DATA = data/observed_stats_extra_line.txt

TESTS = \
	check_array_utils \
	check_stats_utils \
	check_parsing

check_PROGRAMS = \
	check_array_utils \
	check_stats_utils \
	check_parsing

check_array_utils_SOURCES = \
	check_array_utils.c \
	test_utils.c \
	test_utils.h

check_array_utils_CFLAGS = @CHECK_CFLAGS@
check_array_utils_LDADD =  @CHECK_LIBS@ @MATHLDADD@

check_stats_utils_SOURCES = \
	check_stats_utils.c \
	$(top_builddir)/src/array_utils.c \
	test_utils.c \
	test_utils.h

check_stats_utils_CFLAGS = @CHECK_CFLAGS@
check_stats_utils_LDADD =  @CHECK_LIBS@ @MATHLDADD@

check_parsing_SOURCES = \
	check_parsing.c \
	$(top_builddir)/src/array_utils.c \
	$(top_builddir)/src/stats_utils.c \
	test_utils.c \
	test_utils.h

check_parsing_CFLAGS = @CHECK_CFLAGS@
check_parsing_LDADD =  @CHECK_LIBS@
