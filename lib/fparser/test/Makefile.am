bin_PROGRAMS = testfp  \
               testfp1  \
               testfp2  \
               testfp3  \
               testfp4  

testfp_SOURCES = testfp.f90  
testfp1_SOURCES = testfp1.f90  
testfp2_SOURCES = testfp2.f90  
testfp3_SOURCES = testfp3.f90  
testfp4_SOURCES = testfp4.f90  


LDADD = $(top_srcdir)/src/libfparser.a
AM_LDFLAGS = -L$(top_srcdir)/src/ 
AM_FCFLAGS = -I$(top_srcdir)/src/

DEBUG_FLAGS = -g -O0 -p
PROF_FLAGS = -p 

if COMPILER_IS_GFORTRAN
DEBUG_FLAGS += -fbacktrace
DEBUG_FLAGS += -fcheck=all
DEBUG_FLAGS += -Wall
DEBUG_FLAGS += -fbounds-check
DEBUG_FLAGS += -fcheck-array-temporaries 
DEBUG_FLAGS += -pedantic-errors
# F90DEBUG += -std=f2003 -pedantic 

# Some times needed for valgrind
# F90DEBUG += -static 

# To warn about use modules without only keyword
# F90DEBUG += -Wuse-without-only 
endif
    
if COMPILER_IS_IFORT
DEBUG_FLAGS += -check all 
DEBUG_FLAGS += -traceback 
DEBUG_FLAGS += -warn all  
DEBUG_FLAGS += -debug all 
#F90DEBUG += "-Mprof=func"
#F90DEBUG += "-Mprof=line"
endif

if DEBUG
  ENABLE_DEBUG=$(DEBUG_FLAGS)
endif

AM_FCFLAGS += $(ENABLE_DEBUG)



.PHONY: debug

debug:
	make ENABLE_DEBUG="$(DEBUG_FLAGS)"

all: testfp testfp1 testfp2 testfp3 testfp4
