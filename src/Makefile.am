bin_PROGRAMS = gems

gems_SOURCES = Algebra.F90            \
               Fields.f90             \
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
               Pairs.F90                \
	       Bias.f90	              \
               Main.F90               \
               Metadynamics.f90       \
               MPI.F90                \
               NEB.f90                \
               Neighbor.F90           \
               Cells.F90              \
               Output.F90             \
               Programs.f90           \
               Program_Types.F90      \
               Groups.F90	      \
               Graphs.F90	      \
               Calc.F90	              \
               Quasi_Newton.f90       \
               Select_Create.F90      \
               Set_Properties.f90     \
               Strings.f90            \
               Tables.f90             \
               TB.F90                 \
               Variables.F90

    
# zeros...
gems_LDADD =
AM_FCFLAGS =
 
# nodist_gems_SOURCES = $(gems_common_SOURCES)
# gems_common_SOURCES = sprng_f.h
# CLEANFILES = $(gems_common_SOURCES)
#
# $(gems_common_SOURCES):
# 	    cp $(top_builddir)/lib/sprng5/include/$@ $@


AM_FCFLAGS += $(OPENMP_CFLAGS)
AM_FCFLAGS += -DDIM3

# Only compile ReplicaExchange if MPI is ensure
if HAVE_MPI
replicaexchange_obj = ReplicaExchange.o
AM_FCFLAGS += -DHAVE_MPI
gems_SOURCES += ReplicaExchange.f90
ReplicaExchange.o: Metadynamics.o Constants.o Errors.o Program_Types.o Programs.o MPI.o Interaction.o Integration.o $(random_obj)
endif
        
# External software
# =================

# Note on local vs external libraries:
# > We recommend that you avoid using -l options in LDADD or prog_LDADD when
#   referring to libraries built by your package. Instead, write the file name of
#   the library explicitly as in the above cpio example. Use -l only to list
#   third-party libraries. If you follow this rule, the default value of
#   prog_DEPENDENCIES will list all your local libraries and omit the other ones.
# [source](https://www.gnu.org/software/automake/manual/html_node/Linking.html)

# Fortran Preprocesor Templates
# -----------------------------
AM_FCFLAGS += -I$(top_srcdir)/lib/fpt/include/
  
# fparser
# -------
AM_FCFLAGS += -I$(top_srcdir)/lib/fparser/src/
gems_LDADD += $(top_srcdir)/lib/fparser/src/libfparser.a

# sprng
# -----	
if SPRNG
random_obj = Random_sprng.o
gems_SOURCES += Random_sprng.F90 
            
AM_FCFLAGS += -I$(top_srcdir)/lib/sprng5/include 
gems_LDADD += $(top_srcdir)/lib/sprng5/lib/libsprng.a -lstdc++ -lm    
                 
Random_sprng.o: Constants.o  Errors.o

else

random_obj = Random.o
gems_SOURCES += Random.F90 
            
Random.o: Constants.o  Errors.o

endif
        
 
# Debug
# =====
   
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
Main.o: 
Constants.o: 
Strings.o:  Constants.o
Algebra.o:  Constants.o $(random_obj)
Errors.o:   Constants.o Strings.o
Elements.o: Constants.o Strings.o Errors.o

Input_Parsing.o: Constants.o Elements.o Algebra.o 
Program_Types.o: Constants.o Errors.o Elements.o
Groups.o:        Constants.o Elements.o Errors.o Program_Types.o Algebra.o 
Variables.o:     Constants.o Strings.o Errors.o Groups.o 
Cells.o:         Constants.o Algebra.o Program_Types.o Variables.o
Neighbor.o:      Constants.o Algebra.o Program_Types.o Variables.o Cells.o

DDDA.o:   Errors.o   Program_Types.o 
Tables.o: Input_Parsing.o   Program_Types.o DDDA.o

Inq_Properties.o: Constants.o Strings.o Algebra.o    Input_Parsing.o Program_Types.o Neighbor.o    Tables.o
Set_Properties.o: Constants.o Strings.o Algebra.o    Program_Types.o    Inq_Properties.o  
CVS.o: 		  Constants.o Errors.o Algebra.o     Program_Types.o    Inq_Properties.o  
Select_Create.o:  Program_Types.o    Inq_Properties.o Set_Properties.o 

Output.o: Constants.o Errors.o   Input_Parsing.o Program_Types.o Variables.o    DDDA.o   Inq_Properties.o Set_Properties.o     

Bias.o:       Program_Types.o Neighbor.o 
Pairs.o:      Algebra.o    Program_Types.o Neighbor.o Inq_Properties.o Set_Properties.o
Fields.o:     Constants.o Algebra.o   Program_Types.o Neighbor.o   Output.o 
TB.o:         Algebra.o Program_Types.o Neighbor.o Inq_Properties.o Set_Properties.o
Graphs.o:     Constants.o Strings.o   Groups.o Input_Parsing.o     Neighbor.o 
Calc.o:       Constants.o Strings.o   Groups.o Input_Parsing.o     Cells.o  
ForceField.o: Constants.o Algebra.o   Input_Parsing.o Program_Types.o    Output.o
BO.o:         Algebra.o    Program_Types.o Neighbor.o    Inq_Properties.o Set_Properties.o

Bond_Boost.o:   Algebra.o   Program_Types.o   Neighbor.o   Inq_Properties.o Set_Properties.o
Quasi_Newton.o: Algebra.o   Program_Types.o   Inq_Properties.o Set_Properties.o   Output.o 

Integration.o: Constants.o Algebra.o    Program_Types.o  Variables.o Neighbor.o   Inq_Properties.o Set_Properties.o   Bias.o 
Interaction.o: Constants.o    Program_Types.o Neighbor.o    CVS.o   Output.o   TB.o ForceField.o Bias.o Pairs.o  Graphs.o
Checkpoint.o:  Groups.o    Input_Parsing.o    Program_Types.o   Output.o   Integration.o Interaction.o

NEB.o:         Constants.o Algebra.o Errors.o    Program_Types.o   Inq_Properties.o Set_Properties.o   Output.o  Interaction.o Integration.o    

MPI.o:        Constants.o Strings.o
CLInterpreter.o: $(replicaexchange_obj) Hyperdinamics.o Programs.o Metadynamics.o MPI.o Checkpoint.o  Errors.o   Input_Parsing.o  Inq_Properties.o   Integration.o  Interaction.o NEB.o  Neighbor.o  Output.o Program_Types.o  Quasi_Newton.o $(random_obj) Select_Create.o  Set_Properties.o CVS.o  Variables.o Graphs.o Calc.o 
Hyperdinamics.o: Bias.o Integration.o  Quasi_Newton.o  Programs.o  Set_Properties.o     Inq_Properties.o     Program_Types.o	  
Metadynamics.o:  Interaction.o Integration.o  Programs.o   Inq_Properties.o     Program_Types.o	  
Programs.o:      Interaction.o Integration.o Input_Parsing.o Constants.o Output.o Neighbor.o Set_Properties.o Inq_Properties.o  Program_Types.o  Errors.o Checkpoint.o  
        
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
         

clean-local:
	rm -f *.mod

