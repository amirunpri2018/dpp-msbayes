SUBDIRS = src . test
EXTRA_DIST = bootstrap.sh configure.ac

check-local:
	mkdir -p "$(abs_top_builddir)/test/data"
	cp -r "$(abs_top_srcdir)/test/data" "$(abs_top_builddir)/test"

clean-local:
	rm -rf "$(abs_top_builddir)/test/data"
