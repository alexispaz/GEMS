#ifndef _lcg_h
#define _lcg_h

extern "C" {
class LCG : public Sprng
{
 public:

  LCG();
  int init_rng(int, int, int, int);
  ~LCG();
  LCG (const LCG &);
  LCG & operator= (const LCG &);

  int get_rn_int ();
  float get_rn_flt ();
  double get_rn_dbl ();
  int spawn_rng (int nspawned, Sprng ***newgens);
  int get_seed_rng ();
  int free_rng ();
  int pack_rng (char **buffer);
  int unpack_rng (char *packed);
  int print_rng ();

 private:

#ifdef LONG64
  int rng_type;
  unsigned LONG64 seed;
  int init_seed;
  int prime;
  int prime_position;
  int prime_next;
  char *gentype;
  int parameter;
  unsigned LONG64 multiplier;

  inline void multiply();

#else
  int rng_type;
  int seed[2];
  int init_seed;
  int prime;
  int prime_position;
  int prime_next;
  char *gentype;
  int parameter;
  int *multiplier;

  inline void multiply(int *, int *, int *);
  
#endif

  void advance_seed();
};
}

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
