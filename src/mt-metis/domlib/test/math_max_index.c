#include "test.h"

#define N 1000
#define M 197

sint_t test(void) 
{
  sint_t i,j;
  sint_t * ptr = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i%M;
  }
  
  j = sint_max_value(ptr,N);

  TESTEQUALS(j,M-1,PF_SINT_T);
  return 0;
}
