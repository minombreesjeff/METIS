#include "test.h"

#define N 1000

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);
  sint_t * sum = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = sint_rand(0,100);
  }

  sum[0] = ptr[0];
  for (i=1;i<N;++i) {
    sum[i] = ptr[i] + sum[i-1];
  }

  sint_prefixsum_inc(ptr,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(sum[i],ptr[i],PF_SINT_T);
  }
  return 0;
}
