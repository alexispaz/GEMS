#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(sprng,5)
AC_CONFIG_SRCDIR(SRC/lfg/lfg.cpp)
AC_CANONICAL_SYSTEM

# Avoid make to rebuild configure
AM_MAINTAINER_MODE

AC_CONFIG_MACRO_DIRS([m4])

# AM_INIT_AUTOMAKE([subdir-objects])
AM_INIT_AUTOMAKE()

AM_CONFIG_HEADER(config.h:config.in)

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])],
                  with_mpi=$withval, with_mpi=no)

# Checks for programs.
  
AC_PROG_CXX

  AC_PROG_CC
  AC_PROG_CPP

AX_MPI([
  if test x"$with_mpi" != "xno"; then
    use_mpi="yes"
  fi
],[
  use_mpi="no"
  if test x"$with_mpi" = "xyes"; then
    AC_MSG_FAILURE([MPI compiler requested, but could not use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, will not use MPI.])
  fi
])

if test x"$use_mpi" = "xyes"; then
  FC="$MPIFC"
  F77="$MPIF77"
fi   
 

AC_F77_LIBRARY_LDFLAGS 

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

AC_PROG_RANLIB
AC_PATH_PROGS(BASH, bash)

AC_SUBST(sprng_cxxflags)
AC_SUBST(sprng_fflags)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stddef.h stdlib.h string.h sys/time.h])

AC_CHECK_SIZEOF([long int])

CXXFLAGS="$CXXFLAGS -Wno-deprecated"

if test "$ac_cv_sizeof_long_int" = "8"; then
  CXXFLAGS="$CXXFLAGS -DLONG64=long"
fi
			  
AM_CONDITIONAL([USE_MPI], [test $use_mpi = y])

if test $use_mpi = y; then
  MPI_DEF="-DSPRNG_MPI"
else
  MPI_DEF=""
fi

AC_SUBST(MPI_DEF)

AC_ARG_WITH(testu01, [AC_HELP_STRING([--with-testu01], [path to testU01])], 
                  with_testu01=$withval, with_testu01=no)
if test "$with_testu01" = "no"; then
  use_testu01=n
else
  use_testu01=y
  TESTU01_LIBDIR=$with_testu01"/lib"
  TESTU01_INCDIR=$with_testu01"/include"
fi                 

AC_ARG_WITH(testu01lib, [AC_HELP_STRING([--with-testu01lib], [path to testU01 libraries])], 
                  with_testu01lib=$withval, with_testu01lib=no)
AC_ARG_WITH(testu01inc, [AC_HELP_STRING([--with-testu01inc], [path to testU01 include directory])], 
                  with_testu01inc=$withval, with_testu01inc=no)
if test "$with_testu01lib" != "no"; then
  TESTU01_LIBDIR=$with_testu01lib
  if test "$with_testu01inc" != "no"; then
    use_testu01=y
    TESTU01_INCDIR=$with_testu01inc
  fi                 
fi

AC_SUBST(TESTU01_LIBDIR)
AC_SUBST(TESTU01_INCDIR)

AM_CONDITIONAL(USE_TESTU01, test "$use_testu01" = "y")

AC_ARG_WITH(fortran, [AC_HELP_STRING([--with-fortran],[enable fortran compilation])],
                  with_fortran=$withval, with_fortran=yes)                  
                  
if test "$with_fortran" = "yes"; then  
  if test -z "$F77"; then
    use_fortran=y
  fi
else
  use_fortran=n 
fi


AM_CONDITIONAL(compile_FORTRAN, test "$use_fortran" = "y")

AC_CHECK_HEADER(fft.h,fft="yes",fft="no")
if test $fft="yes"; then
  AC_CHECK_FUNC(dzfft2dui, [fft="yes"], [fft="no"])
fi
AM_CONDITIONAL(FFT, test $fft = "yes")

AC_C_BIGENDIAN

dnl deal with Fortran pointer size
AC_CHECK_SIZEOF(int *)
if test "$ac_cv_sizeof_int_p" = "8"; then
  FFPSIZE="-DPOINTER_SIZE=8 -DLONG64=long -DINTEGER_STAR_8"
else
  FFPSIZE=""  
fi

FMANGLE=""

AC_SUBST(FFPSIZE)
AC_SUBST(FMANGLE)

AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(float)

if test "$with_fortran" != "yes"; then  
  F77="echo"
  FORTRANDIR=
else
  FORTRANDIR="F77"

AC_F77_WRAPPERS dnl provides F77_FUNC and  F77_FUNC_ for name mangling
fi

AC_SUBST(FORTRANDIR)

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_CHECK_FUNCS([memset pow sqrt])

AC_PROG_INSTALL
AC_OUTPUT([Makefile 
	     include/Makefile
	     lib/Makefile
	     SRC/Makefile
	     SRC/cmrg/Makefile
	     SRC/lcg/Makefile
	     SRC/lcg64/Makefile
	     SRC/lfg/Makefile
	     SRC/mlfg/Makefile
	     SRC/pmlcg/Makefile
		 ]) 
# AC_OUTPUT([Makefile 
# 	     include/Makefile
# 	     lib/Makefile
# 	     SRC/Makefile
# 	     SRC/cmrg/Makefile
# 	     SRC/lcg/Makefile
# 	     SRC/lcg64/Makefile
# 	     SRC/lfg/Makefile
# 	     SRC/mlfg/Makefile
# 	     SRC/pmlcg/Makefile
# 		 check/Makefile 
# 		 check/lcg/Makefile 
# 		 check/lcg/F77/Makefile 
# 		 check/lfg/Makefile
# 		 check/lfg/F77/Makefile
# 		 check/mlfg/Makefile
# 		 check/mlfg/F77/Makefile
# 		 check/cmrg/Makefile
# 		 check/cmrg/F77/Makefile
# 		 check/lcg64/Makefile
# 		 check/lcg64/F77/Makefile
# 		 check/pmlcg/Makefile
# 		 check/pmlcg/F77/Makefile
# 		 TESTS/Makefile
# 		 TESTS/mpitests/Makefile
#      TESTS/testu01/Makefile
#      EXAMPLES/Makefile
# 		 EXAMPLES/mpisprng/Makefile
# 		 EXAMPLES/F77/Makefile])
