#ifndef _fwrap_h
#define _fwrap_h

/************************************************************************/
/************************************************************************/
/*      Inter-language Naming Convention Problem Solution               */
/*                                                                      */
/*      Note that with different compilers you may find that            */
/*      the linker fails to find certain modules due to the naming      */
/*      conventions implicit in particular compilers. This solution     */
/*      uses autoconf macro AC_F77_WRAPPERS, which defines              */
/*      F77_FUNC_(name, NAME). Turns name into something callable       */
/*      from FORTRAN, typically by appending one or more underscores.   */
/*      The functions are called in Fortran using the first argument    */
/*      of F77_FUNC (which are the same as the defined macro, except    */
/*      without _F77).                                                  */
/*                                                                      */
/************************************************************************/
/************************************************************************/



/* Tests for HAVE_CONFIG_H, which means the configure script has been run properly */
#ifdef HAVE_CONFIG_H
/* Including this to use the macros defined during configure script */
#include "../config.h"
#endif
#ifdef F77_FUNC_
  #define ffree_rng_F77 F77_FUNC_(ffree_rng,FFREE_RNG)
  #define fmake_new_seed_F77 F77_FUNC_(fmake_new_seed,FMAKE_NEW_SEED)
  #define fseed_mpi_F77 F77_FUNC_(fseed_mpi,FSEED_MPI)
  #define finit_rng_F77 F77_FUNC_(finit_rng,FINIT_RNG)
  #define fspawn_rng_F77 F77_FUNC_(fspawn_rng,FSPAWN_RNG)
  #define fget_rn_int_F77 F77_FUNC_(fget_rn_int,FGET_RN_INT)
  #define fget_rn_flt_F77 F77_FUNC_(fget_rn_flt,FGET_RN_FLT)
  #define fget_rn_dbl_F77 F77_FUNC_(fget_rn_dbl,FGET_RN_DBL)
  #define fget_seed_rng_F77 F77_FUNC_(fget_seed_rng,FGET_SEED_RNG)
  #define fpack_rng_F77 F77_FUNC_(fpack_rng,FPACK_RNG)
  #define funpack_rng_F77 F77_FUNC_(funpack_rng,FUNPACK_RNG)
  #define fprint_rng_F77 F77_FUNC_(fprint_rng,FPRINT_RNG)

  #define finit_rng_sim_F77 F77_FUNC_(finit_rng_sim,FINIT_RNG_SIM)
  #define fget_rn_int_sim_F77 F77_FUNC_(fget_rn_int_sim,FGET_RN_INT_SIM)
  #define fget_rn_flt_sim_F77 F77_FUNC_(fget_rn_flt_sim,FGET_RN_FLT_SIM)
  #define fget_rn_dbl_sim_F77 F77_FUNC_(fget_rn_dbl_sim,FGET_RN_DBL_SIM)
  #define finit_rng_simmpi_F77 F77_FUNC_(finit_rng_simmpi,FINIT_RNG_SIMMPI)
  #define fget_rn_int_simmpi_F77 F77_FUNC_(fget_rn_int_simmpi,FGET_RN_INT_SIMMPI)
  #define fget_rn_flt_simmpi_F77 F77_FUNC_(fget_rn_flt_simmpi,FGET_RN_FLT_SIMMPI)
  #define fget_rn_dbl_simmpi_F77 F77_FUNC_(fget_rn_dbl_simmpi,FGET_RN_DBL_SIMMPI)
  #define fpack_rng_simple_F77 F77_FUNC_(fpack_rng_simple,FPACK_RNG_SIMPLE)
  #define funpack_rng_simple_F77 F77_FUNC_(funpack_rng_simple,FUNPACK_RNG_SIMPLE)
  #define fprint_rng_simple_F77 F77_FUNC_(fprint_rng_simple,FPRINT_RNG_SIMPLE)

  #define fcpu_t_F77 F77_FUNC_(fcpu_t,FCPU_T)
#endif

#endif

/***********************************************************************************
* SPRNG (c) 2014 by Florida State University                                       *
*                                                                                  *
* SPRNG is licensed under a                                                        *
* Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. *
*                                                                                  *
* You should have received a copy of the license along with this                   *
* work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.           *
************************************************************************************/
