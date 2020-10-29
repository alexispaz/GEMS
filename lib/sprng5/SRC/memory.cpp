#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

void *_mymalloc(long size, int line, const char *message)
{
  char *temp;

  if(size == 0)
    return NULL;

  temp = (char *) malloc(size);
  
  if(temp == NULL)
    {
      fprintf(stderr,"\nmemory allocation failure in file: %s at line number: %d\n", message, line);
      return NULL;
    }

  return (void *) temp;
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
