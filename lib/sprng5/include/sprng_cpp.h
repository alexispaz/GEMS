#ifndef _sprng_cpp_h_
#define _sprng_cpp_h_

#if __GNUC__ > 3
 #include <stdlib.h>
#endif

#include "sprng.h"
#include "lfg.h"
#include "lcg.h"
#include "lcg64.h"
#include "cmrg.h"
#include "mlfg.h" 
#include "pmlcg.h"

Sprng * SelectType(int typenum);

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
