bin_PROGRAMS = gems

gems_SOURCES = Algebra.F90            \
               Analitical.f90         \
               Bond_Boost.f90         \
               Checkpoint.F90         \
               CLInterpreter.F90      \
               Constants.F90          \
               CVS.F90                \
               DDDA.F90               \
               Elements.F90           \
               Errors.f90             \
               ForceField.F90         \
               Hyperdinamics.f90      \
               Input_Parsing.F90      \
               Inq_Properties.F90     \
               Integration.F90        \
               Interaction.f90        \
               LJ0.f90                \
	       Bias.f90	              \
               Main.F90               \
               Metadynamics.f90       \
               MPI.F90                \
               NEB.f90                \
               Neighbour.F90          \
               Output.F90             \
               Programs.f90           \
               Program_Types.F90      \
               Atoms.F90	      \
               Groups.F90	      \
               Quasi_Newton.f90       \
               Select_Create.F90      \
               Set_Properties.f90     \
               Strings.f90            \
               Tables.f90             \
               TB.F90                 \
               Tersoff.f90            \
               Variables.F90


# Fortran Preprocesor Templates	       
FPTDDS=$(top_srcdir)/lib/fpt/include/
AM_FCFLAGS = -I$(FPTDDS)
include_HEADERS = $(FPTDDS)/prep.inc \
                  $(FPTDDS)/vector_header.inc \
                  $(FPTDDS)/vector_body.inc   \
                  $(FPTDDS)/list_header.inc   \
                  $(FPTDDS)/list_body.inc     \
                  $(FPTDDS)/dlist_header.inc  \
                  $(FPTDDS)/dlist_body.inc    \
                  $(FPTDDS)/cdlist_header.inc \
                  $(FPTDDS)/cdlist_body.inc   \
                  $(FPTDDS)/arrayofptrs_header.inc

# nodist_gems_SOURCES = $(gems_common_SOURCES)
# gems_common_SOURCES = sprng_f.h
# CLEANFILES = $(gems_common_SOURCES)
#
# $(gems_common_SOURCES):
# 	    cp $(top_builddir)/lib/sprng5/include/$@ $@


AM_FCFLAGS += $(OPENMP_CFLAGS)
gems_LDADD =
gems_LDFLAGS =
                           
AM_FCFLAGS += -I$(top_srcdir)/lib/fparser/src/
gems_LDADD += $(top_srcdir)/lib/fparser/src/libfparser.a


AM_FCFLAGS += -DDIM3
      
# Optional compilation of SPRNG for random number generation
if SPRNG
random_obj = Random_sprng.o
gems_SOURCES += Random_sprng.F90 
            
AM_FCFLAGS += -I$(top_srcdir)/lib/sprng5/include 
gems_LDFLAGS += -L$(top_srcdir)/lib/sprng5/lib     
gems_LDADD += -lsprng -lstdc++ -lm    
                 
Random_sprng.o: Constants.o  Errors.o

else

random_obj = Random.o
gems_SOURCES += Random.F90 
            
Random.o: Constants.o  Errors.o

endif
        

# Only compile ReplicaExchange if MPI is ensure
if HAVE_MPI
replicaexchange_obj = ReplicaExchange.o
AM_FCFLAGS += -DHAVE_MPI
gems_SOURCES += ReplicaExchange.f90
ReplicaExchange.o: Metadynamics.o Constants.o Errors.o Program_Types.o Programs.o MPI.o Interaction.o Integration.o $(random_obj)
endif


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

all: gems	       

# Propios
Main.o:  Neighbour.o
Constants.o:     
Strings.o:       Constants.o
Errors.o:        Constants.o Strings.o
Variables.o:	 Errors.o Program_Types.o Constants.o Strings.o
Algebra.o:	 Constants.o $(random_obj)
Elements.o:      Strings.o Constants.o Errors.o
Input_Parsing.o: Elements.o Algebra.o Constants.o
Atoms.o: Elements.o Constants.o 
Groups.o: Constants.o 
Program_Types.o: Elements.o Constants.o Errors.o Input_Parsing.o Atoms.o Groups.o

ForceField.o:    Algebra.o Output.o Program_Types.o Constants.o Input_Parsing.o    
DDDA.o:                 Program_Types.o	 Errors.o
Tables.o:		Program_Types.o  Input_Parsing.o DDDA.o
Inq_Properties.o:   Algebra.o    Constants.o     Strings.o Constants.o    Input_Parsing.o    Program_Types.o  Tables.o
Set_Properties.o:       Inq_Properties.o Algebra.o   Constants.o     Strings.o Constants.o    Program_Types.o
Select_Create.o:        Set_Properties.o    Inq_Properties.o    Program_Types.o  
Quasi_Newton.o:         Set_Properties.o    Inq_Properties.o    Algebra.o     Program_Types.o  Output.o 
Bias.o:		Neighbour.o Program_Types.o  Interaction.o
LJ0.o:		Neighbour.o Set_Properties.o    Inq_Properties.o     Algebra.o    Program_Types.o
TB.o:		Neighbour.o Set_Properties.o    Inq_Properties.o     Algebra.o    Program_Types.o
Tersoff.o:      Neighbour.o Set_Properties.o    Inq_Properties.o     Algebra.o    Program_Types.o
Bond_Boost.o:   Neighbour.o Set_Properties.o    Inq_Properties.o     Algebra.o    Program_Types.o
BO.o:		Neighbour.o Set_Properties.o    Inq_Properties.o     Algebra.o    Program_Types.o
Analitical.o:   Output.o Neighbour.o Constants.o Algebra.o    Program_Types.o
Integration.o:	 Variables.o Bias.o Algebra.o     Neighbour.o Set_Properties.o    Inq_Properties.o     Constants.o   Program_Types.o
Neighbour.o:   Variables.o   Algebra.o     Set_Properties.o   Program_Types.o Constants.o Inq_Properties.o
Interaction.o: Neighbour.o TB.o ForceField.o  Constants.o   Program_Types.o Output.o Neighbour.o CVS.o
Output.o:  	        DDDA.o     Set_Properties.o     Inq_Properties.o    Constants.o    Program_Types.o Errors.o Input_Parsing.o 
NEB.o:         		Interaction.o Integration.o  Set_Properties.o     Inq_Properties.o     Program_Types.o   Constants.o Errors.o Output.o Algebra.o  

Checkpoint.o:    Input_Parsing.o Program_Types.o  Integration.o Interaction.o Output.o


MPI.o: 		 Constants.o Strings.o
CLInterpreter.o: $(replicaexchange_obj) Hyperdinamics.o Programs.o Metadynamics.o MPI.o Analitical.o  Checkpoint.o  Errors.o  ForceField.o  Input_Parsing.o  Inq_Properties.o   Integration.o  Interaction.o  Bias.o LJ0.o  NEB.o  Neighbour.o  Output.o Program_Types.o  Quasi_Newton.o $(random_obj) Select_Create.o  Set_Properties.o TB.o  Tersoff.o  CVS.o  Variables.o
Hyperdinamics.o: Bias.o Integration.o  Quasi_Newton.o  Programs.o  Set_Properties.o     Inq_Properties.o     Program_Types.o	  
Metadynamics.o:  Interaction.o Integration.o  Programs.o   Inq_Properties.o     Program_Types.o	  
Programs.o:      Interaction.o Integration.o Input_Parsing.o Constants.o Output.o Neighbour.o Set_Properties.o Inq_Properties.o  Program_Types.o  Errors.o Checkpoint.o 
        
CVS.o: Algebra.o  Inq_Properties.o  Program_Types.o Constants.o  Errors.o

if COMPILER_IS_IFORT
  
# fast is the equivalent of -ipo, -O3, -no-prec-div, -static, and -xHost
# F90FAST += -fast

# Una forma facil de decidir cual poner aca es grep flags /proc/cpuinfo y ver
# la mas avanzada

#FCFLAGS += -xHost       # La mejor del Host actual
#F90FAST += -xCORE-AVX2  # a future Intel 
#F90FAST += -xCORE-AVX-I # a future Intel 
#F90FAST += -xAVX        # 2nd Generation Intel Core 
#F90FAST += -xSSE4.2     # Intel Core i7,i5,i3   Intel Xeon 55XX,56XX,75XX 
#F90FAST += -xSSE4.1     # Intel Xeon 74XX  
#F90FAST += -xSSE4.1     # Quad-Core Intel Xeon 54XX, 33XX  
#F90FAST += -xSSE4.1     # Dual-Core Intel Xeon 52XX, 31XX 
#F90FAST += -xSSE4.1     # Intel Core 2 Extreme 9XXX 
#F90FAST += -xSSE4.1     # Intel Core 2 Quad 9XXX 
#F90FAST += -xSSE4.1     # Intel Core 2 Duo 8XXX, E7200
#F90FAST += -xSSSE3      # Quad-Core Intel Xeon 73XX, 53XX, 32XX 
#F90FAST += -xSSSE3      # Dual-Core Intel Xeon 72XX, 53XX, 51XX, 30XX 
#F90FAST += -xSSSE3      # Intel Core 2 Extreme 7XXX, 6XXX 
#F90FAST += -xSSSE3      # Intel Core 2 Quad 6XXX 
#F90FAST += -xSSSE3      # Intel Core 2 Duo 7XXX (except E7200), 6XXX, 5XXX, 4XXX 
#F90FAST += -xSSSE3      # Intel Core 2 Solo 2XXX 
#F90FAST += -xSSSE3      # Intel Pentium dual-core  E2XXX, T23XX 
#F90FAST += -xSSE3_ATOM  # Intel Atom 
#F90FAST += -xSSE3       # Dual-Core Intel Xeon 70XX, 71XX, 50XX, (ULV and LV) 1.66, 2.0, 2.162.8
#F90FAST += -xSSE3       # Intel Xeo  with SSE3 instruction set support
#F90FAST += -xSSE3       # Intel Core Duo, Solo
#F90FAST += -xSSE3       # Intel Pentium dual-core  T21XX, T20XX 
#F90FAST += -xSSE3       # Intel Pentium  Extreme Edition, D
#F90FAST += -xSSE3       # Intel Pentium 4  with SSE3 instruction set support
#F90FAST += -xSSE2       # Intel Xeon 
#F90FAST += -xSSE2       # Intel Pentium 4, M
#F90FAST += -xSSE2       # IA32    Intel Pentium III, II, I 
        

# ifdef DMKL 
# LIBPATHS += -L/opt/intel/composerxe/mkl/lib/intel64/
# INCLUDES += -I/opt/intel/composerxe/mkl/include/intel64/lp64/
# LIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64
# LIBS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
# endif
           


endif




if COMPILER_IS_GFORTRAN
# No limita el tamaño de las lineas
# FCFLAGS += -ffree-line-length-none

#
# ifdef DMKL
# LIBPATHS += -L/usr/lib64
# LIBS += -llapack
# LIBS += -lblas -lcblas
# # Threaded version may work as well
# # LIB_BLAS = -lptf77blas -lptcblas -latlas
# endif
#       
endif
         



