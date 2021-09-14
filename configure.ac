#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([gems],[m4_esyscmd_s([git describe --always --abbrev=4 --dirty])],[alexis.paz@gmail.com])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR(src/)

# AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_INIT_AUTOMAKE([subdir-objects foreign])

# Avoid make to rebuild configure
AM_MAINTAINER_MODE

# Set OpenMP flags
# : ${enable_openmp=no}
AC_OPENMP

# Set AR
AM_PROG_AR

# Switch for MPI
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found, MPI is not used. Default: no])
],,[with_mpi=no])
 
# Set FCFLAGS and FC
: ${FCFLAGS=""} # Override default (i.e. "-g -O2")
# AC_PROG_FC # Replaced by AX_PROG_FC_MPI that call it if needed
AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but could not use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, will not use MPI.])
  fi
])

# # Set MPILIBS and MPICC, MPICXX, MPIF77, or MPIFC (depending upon the current language)
# AX_MPI([
#   if test x"$with_mpi" != "xno"; then
#     use_mpi="yes"
#   fi
# ],[
#   use_mpi="no"
#   if test x"$with_mpi" = "xyes"; then
#     AC_MSG_FAILURE([MPI compiler requested, but could not use MPI.])
#   else
#     AC_MSG_WARN([No MPI compiler found, will not use MPI.])
#   fi
# ])

    
# Debug options
AC_ARG_ENABLE([debug],
  AS_HELP_STRING([--enable-debug],[Turn on debugging]),
  [case "${enableval}" in
    yes) debug=true ;;
    no)  debug=false ;;
    *) AC_MSG_ERROR([bad value "${enableval}" for --enable-debug]) ;;
  esac],enable_debug=no)


# External software
AC_ARG_WITH([sprng], 
  AS_HELP_STRING([--with-sprng],[Use sprng library for random number generation]),,with_sprng=no)
AC_CONFIG_SUBDIRS([
    lib/sprng5
    lib/fparser
])
 
# Usefull variables to use inside the Makefile.am files
AM_CONDITIONAL([COMPILER_IS_IFORT], [test x"$FC" = xifort])
AM_CONDITIONAL([COMPILER_IS_GFORTRAN], [test x"$FC" = xgfortran])
AM_CONDITIONAL([HAVE_MPI], [test x"$use_mpi" = "xyes"])
AM_CONDITIONAL([DEBUG], [test "x$enable_debug" = "xyes"])
AM_CONDITIONAL([SPRNG], [test "x$with_sprng" = "xyes"])

# To check for a compile flag
# AX_CHECK_COMPILE_FLAG([-std=c11],
#   [AX_APPEND_FLAG([-std=c11])],
#   [AC_MSG_WARN([-std=c11 not supported, you may have to set CFLAGS to enable C11 support.])
# ])

AC_CONFIG_FILES([
 Makefile
 lib/fparser/Makefile
 src/Makefile
])


# A precompiler variable
# AC_DEFINE([MYDEFINE],[SOMEVALUE],[Description: Make #define MYDEFINE SOMEVALUE])

AC_OUTPUT  
