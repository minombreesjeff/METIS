#include "test.h"

#define N 1000

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);
  sint_t * b = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i;
  }

  sint_copy(b,ptr,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(b[i],i,PF_SINT_T);
  }
  return 0;
}
