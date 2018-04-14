#include "test.h"

#define N 1000
#define M 197

sint_t test(void) 
{
  sint_t i,j,k;
  sint_t * ptr = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i%M;
  }
  
  j = sint_sum(ptr,N);

  k = 0;
  for (i=0;i<N;++i) {
    k += ptr[i]; 
  }

  TESTEQUALS(j,k,PF_SINT_T);
  return 0;
}
