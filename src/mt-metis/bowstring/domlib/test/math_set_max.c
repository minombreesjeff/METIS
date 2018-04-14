#include "test.h"

#define N 1000

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i;
  }

  sint_set_max(ptr,N-11,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],i-10,PF_SINT_T);
  }
  return 0;
}
