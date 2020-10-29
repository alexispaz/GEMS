#define MAXVAL "152D02C7E14AF6800000" /*maximum # of independent streams*/

struct BIGNUM_ARRAY_TYPE 
{
  long size;
  BigNum *list;
};

#define  OP_SIZE    2
#define  RNGBITS       3
#define  SHIFT      29
#define  MASK       0X1FFFFFFF
#define  MAGIC_NUM  "222222222222222"
#define  MAGIC_DEN  "60454F554C0000"
#define  PRIM       37
#define  POWER_N    61

/* Values pertain to this particular parameter: 2^61-1 as modulus*/
struct BIGNUM_ARRAY_TYPE init_factors()
{
  struct BIGNUM_ARRAY_TYPE factors;

  factors.size = 12;
  factors.list = new BigNum[12];
	
  factors.list[0] = Set_ui(2ul);
  factors.list[1] = Set_ui(3ul);
  factors.list[2] = Set_ui(5ul);
  factors.list[3] = Set_ui(7ul);
  factors.list[4] = Set_ui(11ul);
  factors.list[5] = Set_ui(13ul);
  factors.list[6] = Set_ui(31ul);
  factors.list[7] = Set_ui(41ul);
  factors.list[8] = Set_ui(61ul);
  factors.list[9] = Set_ui(151ul);
  factors.list[10] = Set_ui(331ul);
  factors.list[11] = Set_ui(1321ul);

  return (factors);
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
