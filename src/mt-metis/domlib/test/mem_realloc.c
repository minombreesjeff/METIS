#include "test.h"

#define N 1000
#define M 700

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i;
  }

  ptr = sint_realloc(ptr,M);

  for (i=0;i<M;++i) {
    TESTEQUALS(ptr[i],i,PF_SINT_T);
  }
  return 0;
}
