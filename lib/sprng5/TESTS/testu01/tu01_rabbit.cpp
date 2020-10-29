 #include <iostream>
 #include <cstdio> 
 #include <cstdlib>
 
 #include "sprng_cpp.h"
 extern "C" {
 #include "unif01.h"
 #include "bbattery.h"
 }

 using namespace std;
 
 Sprng * ptr;
 
 /* Returns generated double-precision numbers for testU01 generator */
 double genRanDbl(void) {
   return ptr->sprng();
 }

 
 int main(int argc, char* argv[]) {
  int streamnum, nstreams;
  int gtype, seed, param;
  int nbits;  
  unif01_Gen *gen;

  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " genType seed[optional] parameter[optional] nbits[optional, default=500]\n";
    return 1;
  }
  
  gtype = atoi(argv[1]);
  
  if (argc > 2) {
    seed = atoi(argv[2]);
  } else {
    seed = 985456376;
  }
  
  if (argc > 3) {
    param = atoi(argv[3]);
  } else {
    param = SPRNG_DEFAULT;
  }
  
  if (argc > 4) {
    nbits = atoi(argv[4]);
  } else {
    nbits = 500;
  }
  
  streamnum = 0;
  nstreams = 1;
 
  ptr = SelectType(gtype);
  ptr->init_sprng(streamnum, nstreams, seed, param);

  printf("\n --------------------------------------------------------\n");
  printf(" Printing information about stream:\n");
  ptr->print_sprng();
 
  gen = unif01_CreateExternGen01 ("SPRNG Double", genRanDbl); /* Create testU01 generator */
  
  bbattery_Rabbit(gen, nbits);
  
  unif01_DeleteExternGen01 (gen); /* Delete testU01 generator */
  ptr->free_sprng();  /* free memory used to store stream state */
 
  return 0;
}
 