#include "sprng_cpp.h"

#include <iostream>

using namespace std;

Sprng * SelectType(int typenum)
{
  switch (typenum)
    {
    case 0: return new LFG;
    case 1: return new LCG;
    case 2: return new LCG64;
    case 3: return new CMRG;
    case 4: return new MLFG;
    case 5: return new PMLCG;
    default: cerr << "Invalid generator type number.\n"; exit(EXIT_FAILURE);
    }
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
