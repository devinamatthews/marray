VPATH += $(srcdir)

AM_CPPFLAGS = -I$(srcdir)/include -I$(srcdir)/external/googletest/include

bin_PROGRAMS = test
test_SOURCES = test/iterator.cpp test/marray.cpp test/utility.cpp
test_LDADD = external/googletest/lib/libgtest_main.la external/googletest/lib/libgtest.la

SUBDIRS = external/googletest
