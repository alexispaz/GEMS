SUBDIRS = .

LIBDIR=../../lib
SRCDIR=../../SRC
INCLUDEDIR=../../include
EXTRA_DIST = fft.cpp

GEN=sprng
CLEANFILES=*.$(GEN)

bin_PROGRAMS = equidist.tmp perm.tmp serial.tmp poker.tmp collisions.tmp gap.tmp maxt.tmp sum.tmp coupon.tmp runs.tmp random_walk.tmp wolff.tmp metropolis.tmp

equidist_tmp_SOURCES = equidist.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
equidist_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
equidist_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
equidist_tmp_LDFLAGS = 

perm_tmp_SOURCES = perm.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
perm_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
perm_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm
perm_tmp_LDFLAGS = 

serial_tmp_SOURCES = serial.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
serial_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
serial_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
serial_tmp_LDFLAGS = 

poker_tmp_SOURCES = poker.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h stirling.cpp $(SRCDIR)/sprng_cpp.h
poker_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
poker_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
poker_tmp_LDFLAGS = 

collisions_tmp_SOURCES = collisions.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
collisions_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
collisions_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm
collisions_tmp_LDFLAGS = 

gap_tmp_SOURCES = gap.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
gap_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
gap_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
gap_tmp_LDFLAGS = 

maxt_tmp_SOURCES = maxt.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
maxt_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
maxt_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
maxt_tmp_LDFLAGS = 

sum_tmp_SOURCES = sum.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
sum_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
sum_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm
sum_tmp_LDFLAGS = 

coupon_tmp_SOURCES = coupon.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h stirling.cpp $(SRCDIR)/sprng_cpp.h
coupon_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
coupon_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
coupon_tmp_LDFLAGS = 

runs_tmp_SOURCES = runs.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
runs_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
runs_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
runs_tmp_LDFLAGS = 

random_walk_tmp_SOURCES = random_walk.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
random_walk_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
random_walk_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
random_walk_tmp_LDFLAGS = 

wolff_tmp_SOURCES = wolff.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
wolff_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
wolff_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
wolff_tmp_LDFLAGS = 

metropolis_tmp_SOURCES = metropolis.cpp chisquare.cpp communicate.cpp init_tests.cpp util.h util.cpp tests.h $(SRCDIR)/sprng_cpp.h
metropolis_tmp_CXXFLAGS = -I$(INCLUDEDIR) -DSPRNG_MPI 
metropolis_tmp_LDADD = -L$(LIBDIR) -l$(GEN) -lm 
metropolis_tmp_LDFLAGS = 

bin_SCRIPTS = sprngtest
sprngtest:
	@rm -f *.o
	@mv equidist.tmp equidist.$(GEN)
	@mv perm.tmp perm.$(GEN)
	@mv serial.tmp serial.$(GEN)
	@mv poker.tmp poker.$(GEN)
	@mv collisions.tmp collisions.$(GEN)
	@mv gap.tmp gap.$(GEN)
	@mv maxt.tmp maxt.$(GEN)
	@mv sum.tmp sum.$(GEN)
	@mv coupon.tmp coupon.$(GEN)
	@mv runs.tmp runs.$(GEN)
	@mv random_walk.tmp random_walk.$(GEN)
	@mv wolff.tmp wolff.$(GEN)
	@mv metropolis.tmp metropolis.$(GEN)

