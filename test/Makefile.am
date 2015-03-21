bin_PROGRAMS = test
test_SOURCES = test.cxx

VPATH += $(srcdir)/..
CXXFLAGS += -std=c++11 -I$(srcdir)/.. -Wno-narrowing
CLEANFILES = test.cxx

test.cxx: marray.hpp utility.hpp iterator.hpp test.cxx.header test.cxx.footer
	cat test.cxx.header > test.cxx
	$(CPP) $(CXXFLAGS) -DMARRAY_TEST\(...\)=MARRAY_TEST_BEGIN\{__VA_ARGS__\}MARRAY_TEST_END $< | $(GREP) -o 'MARRAY_TEST_BEGIN.*MARRAY_TEST_END' | $(SED) 's/MARRAY_TEST_BEGIN//g' | $(SED) 's/MARRAY_TEST_END//g' >> test.cxx
	cat test.cxx.footer >> test.cxx
