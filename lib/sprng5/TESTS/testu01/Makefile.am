LIBDIR=../../lib
SRCDIR=../../SRC
INCLUDEDIR=../../include

GEN=sprng
CLEANFILES=*.$(GEN)

bin_PROGRAMS = alphabit.tmp bigCrush.tmp crush.tmp FIPS_140_2.tmp rabbit.tmp smallCrush.tmp


alphabit_tmp_SOURCES = tu01_alphabit.cpp $(SRCDIR)/sprng_cpp.h
alphabit_tmp_CXXFLAGS = -I$(INCLUDEDIR) -I$(TESTU01_INCDIR)
alphabit_tmp_LDADD = -L$(LIBDIR) -L$(TESTU01_LIBDIR) -l$(GEN) -ltestu01 -lm

bigCrush_tmp_SOURCES = tu01_bigCrush.cpp $(SRCDIR)/sprng_cpp.h
bigCrush_tmp_CXXFLAGS = -I$(INCLUDEDIR) -I$(TESTU01_INCDIR)
bigCrush_tmp_LDADD = -L$(LIBDIR) -L$(TESTU01_LIBDIR) -l$(GEN) -ltestu01 -lm 

crush_tmp_SOURCES = tu01_crush.cpp $(SRCDIR)/sprng_cpp.h
crush_tmp_CXXFLAGS = -I$(INCLUDEDIR) -I$(TESTU01_INCDIR)
crush_tmp_LDADD = -L$(LIBDIR) -L$(TESTU01_LIBDIR) -l$(GEN) -ltestu01 -lm 

FIPS_140_2_tmp_SOURCES = tu01_FIPS_140_2.cpp $(SRCDIR)/sprng_cpp.h
FIPS_140_2_tmp_CXXFLAGS = -I$(INCLUDEDIR) -I$(TESTU01_INCDIR)
FIPS_140_2_tmp_LDADD = -L$(LIBDIR) -L$(TESTU01_LIBDIR) -l$(GEN) -ltestu01 -lm 

rabbit_tmp_SOURCES = tu01_rabbit.cpp $(SRCDIR)/sprng_cpp.h
rabbit_tmp_CXXFLAGS = -I$(INCLUDEDIR) -I$(TESTU01_INCDIR)
rabbit_tmp_LDADD = -L$(LIBDIR) -L$(TESTU01_LIBDIR) -l$(GEN) -ltestu01 -lm 

smallCrush_tmp_SOURCES = tu01_smallCrush.cpp $(SRCDIR)/sprng_cpp.h
smallCrush_tmp_CXXFLAGS = -I$(INCLUDEDIR) -I$(TESTU01_INCDIR)
smallCrush_tmp_LDADD = -L$(LIBDIR) -L$(TESTU01_LIBDIR) -l$(GEN) -ltestu01 -lm 

bin_SCRIPTS = sprngtestu01
sprngtestu01:
	@rm -f *.o
	@cp alphabit.tmp alphabit.$(GEN)
	@cp bigCrush.tmp bigCrush.$(GEN)
	@cp crush.tmp crush.$(GEN)
	@cp FIPS_140_2.tmp FIPS_140_2.$(GEN)
	@cp rabbit.tmp rabbit.$(GEN)
	@cp smallCrush.tmp smallCrush.$(GEN)
	@rm alphabit.tmp
	@rm bigCrush.tmp
	@rm crush.tmp
	@rm FIPS_140_2.tmp
	@rm rabbit.tmp
	@rm smallCrush.tmp
#endif