#ifndef _primes_64_h_
#define _primes_64_h_

int getprime_64 (int need, unsigned int *prime_array, int offset);
 
#define MAXPRIME 3037000501U  /* largest odd # < sqrt(2)*2^31+2 */
#define MINPRIME 55108   /* sqrt(MAXPRIME) */
#define MAXPRIMEOFFSET 146138719U /* Total number of available primes */

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
