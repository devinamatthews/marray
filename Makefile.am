VPATH += $(srcdir)

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(srcdir)/include -I$(srcdir)/external/googletest/include

bin_PROGRAMS = bin/test
bin_test_SOURCES = test/miterator.cxx test/marray.cxx test/utility.cxx
bin_test_LDADD = external/googletest/lib/libgtest_main.la external/googletest/lib/libgtest.la

SUBDIRS = external/googletest
