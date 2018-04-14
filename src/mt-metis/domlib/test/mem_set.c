#include "test.h"

#define N 1000
#define M 0x728fa342

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  sint_set(ptr,M,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],M,PF_SINT_T);
  }
  return 0;
}
