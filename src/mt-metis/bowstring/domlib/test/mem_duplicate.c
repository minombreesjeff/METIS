#include "test.h"

#define N 1000

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i;
  }

  sint_t * dup = sint_duplicate(ptr,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(dup[i],ptr[i],PF_SINT_T);
  }
  return 0;
}
