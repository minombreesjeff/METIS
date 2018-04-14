#include "test.h"

#define N ((sint_t)0x0000F0F0UL)
#define NR ((sint_t)0x0F0F0000UL)
#define M ((sint_t)0x80001000UL)
#define MR ((sint_t)0x00080001UL)

sint_t test(void) 
{
  TESTEQUALS(NR,sint_reversebits(N),PF_SINT_T);
  TESTEQUALS(MR,sint_reversebits(M),PF_SINT_T);
  TESTEQUALS(N,sint_reversebits(NR),PF_SINT_T);
  TESTEQUALS(M,sint_reversebits(MR),PF_SINT_T);
  return 0;
}
