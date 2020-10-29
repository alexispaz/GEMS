/*************************************************************************/
/*************************************************************************/
/* Header file for arithmetic on large rationals                         */
/*                                                                       */
/* Author: J. Ren,                                                       */
/*            Florida State Unversity                                    */
/* E-Mail: ren@csit.fsu.edu                                              */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

#ifndef BIGRAT_H
#define BIGRAT_H

#include "bignum.h"

class BigRat
{
 public:

  BigRat();

  BigNum br_get_num();
  BigNum br_get_den();

  void br_set_num(const BigNum &);
  void br_set_den(const BigNum &);

  BigNum numerator;
  BigNum denominator;
};

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
