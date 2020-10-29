
SUBDIRS = $(FORTRANDIR) .

EXTRA_DIST=mlfg.data .insertlib
INCLUDEDIR = $(top_srcdir)/SRC

LIBDIR=../../lib
CHKDIR=../../check

bin_PROGRAMS = check.cmlfg check_sim.cmlfg time.cmlfg
check_cmlfg_SOURCES = check_gen.cpp 
check_cmlfg_DEPENDENCIES = $(LIBDIR)/libsprng.a
check_cmlfg_CXXFLAGS = -I$(INCLUDEDIR)
check_cmlfg_LDADD = -L$(LIBDIR) -lsprng
check_cmlfg_LDFLAGS = 
check_sim_cmlfg_SOURCES = check_gen_simple.cpp
check_sim_cmlfg_DEPENDENCIES = $(LIBDIR)/libsprng.a
check_sim_cmlfg_CXXFLAGS = -I$(INCLUDEDIR)
check_sim_cmlfg_LDADD = -L$(LIBDIR) -lsprng
check_sim_cmlfg_LDFLAGS = 
time_cmlfg_SOURCES = timing.cpp
time_cmlfg_DEPENDENCIES = $(LIBDIR)/libsprng.a
time_cmlfg_LDADD = -L$(LIBDIR) -lsprng
time_cmlfg_CXXFLAGS = -I$(INCLUDEDIR)
time_cmlfg_LDFLAGS = 

bin_SCRIPTS = checkmlfg timemlfg
CLEANFILES = $(bin_SCRIPTS)

checkmlfg: $(CHKDIR)/insertlib
	rm -f ./checkmlfg
	@$(CHKDIR)/insertlib mlfg ./checkmlfg "Checking Default C++ interface:" "check.cmlfg < ./mlfg.data 2> /dev/null" ./check.cmlfg
	@$(CHKDIR)/insertlib mlfg ./checkmlfg "Checking Simple C++ interface:" "check_sim.cmlfg < ./mlfg.data 2> /dev/null" ./check_sim.cmlfg
	chmod u+x ./checkmlfg

timemlfg: $(CHKDIR)/insertlib
	rm -f ./timemlfg
	@$(CHKDIR)/insertlib mlfg ./timemlfg "Timing C++ interface" ./time.cmlfg ./time.cmlfg
	chmod u+x ./timemlfg
