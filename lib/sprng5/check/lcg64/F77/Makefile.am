
SUBDIRS = .

# EXTRA_DIST=lcg.data .insertlib
INCLUDEDIR = $(top_srcdir)/SRC

LIBDIR=../../../lib
CHKDIR=../../../check

bin_PROGRAMS = check.flcg check_sim.flcg time.flcg
check_flcg_SOURCES = check_genf.F 
check_flcg_DEPENDENCIES = $(LIBDIR)/libsprng.a
check_flcg_FFLAGS = $(FFPSIZE) -I$(INCLUDEDIR)
check_flcg_LDADD = -L$(LIBDIR) -lsprng -lstdc++ 
check_flcg_LDFLAGS = 
check_sim_flcg_SOURCES = check_genf_simple.F 
check_sim_flcg_DEPENDENCIES = $(LIBDIR)/libsprng.a
check_sim_flcg_FFLAGS = $(FFPSIZE) -I$(INCLUDEDIR) 
check_sim_flcg_LDADD = -L$(LIBDIR) -lsprng -lstdc++ 
check_sim_flcg_LDFLAGS = 
time_flcg_SOURCES = timingf.F
time_flcg_DEPENDENCIES = $(LIBDIR)/libsprng.a
time_flcg_FFLAGS = $(FFPSIZE) -I$(INCLUDEDIR)
time_flcg_LDADD = -L$(LIBDIR) -lsprng -lstdc++ 
time_flcg_LDFLAGS = 


bin_SCRIPTS = checklcg timelcg
CLEANFILES = $(bin_SCRIPTS)

checklcg: $(CHKDIR)/insertlib
	rm -f ./checklcg
	@$(CHKDIR)/insertlib lcg ./checklcg "Checking Default FORTRAN interface:" "check.flcg < ./lcg64.data 2> /dev/null" ./check.flcg
	@$(CHKDIR)/insertlib lcg ./checklcg "Checking Simple FORTRAN interface:" "check_sim.flcg < ./lcg64.data 2> /dev/null" ./check_sim.flcg
	chmod u+x ./checklcg

timelcg: $(CHKDIR)/insertlib
	rm -f ./timelcg
	@$(CHKDIR)/insertlib lcg ./timelcg "Timing FORTRAN interface" ./time.flcg ./time.flcg
	chmod u+x ./timelcg
