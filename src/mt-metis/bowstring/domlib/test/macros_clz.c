#include "test.h"

#define N 0x0000000000000001UL
#define N_L 63UL
#define M 0x0UL
#define M_L 64UL
#define O 0x0000F00004000FFFUL
#define O_L 16UL
#define A 0x0000F0F0U
#define A_L 16UL
#define B 0x00000000U
#define B_L 32UL
#define C 0xF0000000U
#define C_L 0UL
#define D 0x805555555555555FUL
#define D_L 0UL

sint_t test(void) 
{
  TESTEQUALS(dl_clz(N),N_L,PF_SIZE_T);
  TESTEQUALS(dl_clz(M),M_L,PF_SIZE_T);
  TESTEQUALS(dl_clz(O),O_L,PF_SIZE_T);
  TESTEQUALS(dl_clz(A),A_L,PF_SIZE_T);
  TESTEQUALS(dl_clz(B),B_L,PF_SIZE_T);
  TESTEQUALS(dl_clz(C),C_L,PF_SIZE_T);
  TESTEQUALS(dl_clz(D),D_L,PF_SIZE_T);
  return 0;
}
