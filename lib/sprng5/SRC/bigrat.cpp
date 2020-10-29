/*************************************************************************/
/*************************************************************************/
/* Implementation file for arithmetic on large rationals                 */
/*                                                                       */
/* Author: J. Ren,                                                       */
/*            Florida State Unversity                                    */
/* E-Mail: ren@csit.fsu.edu                                              */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

#include "bigrat.h"

BigRat::BigRat()
{
  numerator = 0ul;
  denominator = 1ul;
}

BigNum BigRat::br_get_num()
{
  return numerator;
}

BigNum BigRat::br_get_den()
{
  return denominator;
}

void BigRat::br_set_num(const BigNum & num)
{
  numerator = num;
}

void BigRat::br_set_den(const BigNum & den)
{
  denominator = den;
}


/***********************************************************************************
* SPRNG (c) 2014 by Florida State University                                       *
*                                                                                  *
* SPRNG is licensed under a                                                        *
* Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. *
*                                                                                  *
* You should have received a copy of the license along with this                   *
* work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.           *
************************************************************************************/
