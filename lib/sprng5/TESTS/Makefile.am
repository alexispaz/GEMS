if USE_MPI
  SUBD = . mpitests
else
  SUBD = .
endif

if USE_TESTU01
  SUBDIRS = $(SUBD) testu01
else
  SUBDIRS = $(SUBD)
endif

LIBDIR=../lib
SRCDIR=../SRC
INCLUDEDIR=../include
EXTRA_DIST = fft.cpp

GEN=sprng
CLEANFILES=*.$(GEN)

bin_PROGRAMS = equidist.tmp perm.tmp serial.tmp poker.tmp collisions.tmp gap.tmp maxt.tmp sum.tmp coupon.tmp runs.tmp random_walk.tmp wolff.tmp metropolis.tmp wolfftest.tmp wolffind.tmp

equidist_tmp_SOURCES = equidist.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
equidist_tmp_CXXFLAGS = -I$(INCLUDEDIR)
equidist_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

perm_tmp_SOURCES = perm.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
perm_tmp_CXXFLAGS = -I$(INCLUDEDIR)
perm_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm

serial_tmp_SOURCES = serial.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
serial_tmp_CXXFLAGS = -I$(INCLUDEDIR)
serial_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

poker_tmp_SOURCES = poker.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h stirling.cpp $(SRCDIR)/sprng_cpp.h
poker_tmp_CXXFLAGS = -I$(INCLUDEDIR)
poker_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

collisions_tmp_SOURCES = collisions.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
collisions_tmp_CXXFLAGS = -I$(INCLUDEDIR)
collisions_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

gap_tmp_SOURCES = gap.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
gap_tmp_CXXFLAGS = -I$(INCLUDEDIR)
gap_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

maxt_tmp_SOURCES = maxt.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
maxt_tmp_CXXFLAGS = -I$(INCLUDEDIR)
maxt_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

sum_tmp_SOURCES = sum.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
sum_tmp_CXXFLAGS = -I$(INCLUDEDIR)
sum_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm

coupon_tmp_SOURCES = coupon.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h stirling.cpp $(SRCDIR)/sprng_cpp.h
coupon_tmp_CXXFLAGS = -I$(INCLUDEDIR)
coupon_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm

runs_tmp_SOURCES = runs.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
runs_tmp_CXXFLAGS = -I$(INCLUDEDIR)
runs_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

random_walk_tmp_SOURCES = random_walk.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
random_walk_tmp_CXXFLAGS = -I$(INCLUDEDIR)
random_walk_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm

wolff_tmp_SOURCES = wolff.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
wolff_tmp_CXXFLAGS = -I$(INCLUDEDIR)
wolff_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm

wolffind_tmp_SOURCES = wolffind.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
wolffind_tmp_CXXFLAGS = -I$(INCLUDEDIR)
wolffind_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

wolfftest_tmp_SOURCES = wolfftest.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
wolfftest_tmp_CXXFLAGS = -I$(INCLUDEDIR)
wolfftest_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm

metropolis_tmp_SOURCES = metropolis.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
metropolis_tmp_CXXFLAGS = -I$(INCLUDEDIR)
metropolis_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 

bin_SCRIPTS = sprngtest
sprngtest:
	@rm -f *.o
	@cp equidist.tmp equidist.$(GEN)
	@cp perm.tmp perm.$(GEN)
	@cp serial.tmp serial.$(GEN)
	@cp poker.tmp poker.$(GEN)
	@cp collisions.tmp collisions.$(GEN)
	@cp gap.tmp gap.$(GEN)
	@cp maxt.tmp maxt.$(GEN)
	@cp sum.tmp sum.$(GEN)
	@cp coupon.tmp coupon.$(GEN)
	@cp runs.tmp runs.$(GEN)
	@cp random_walk.tmp random_walk.$(GEN)
	@cp wolff.tmp wolff.$(GEN)
	@cp wolffind.tmp wolffind.$(GEN)
	@cp wolfftest.tmp wolfftest.$(GEN)
	@cp metropolis.tmp metropolis.$(GEN)

