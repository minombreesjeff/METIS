#include "test.h"

#define N 1000000
#define M 1971

sint_t test(void) 
{
  sint_t i;
  real_t j,k;
  real_t * ptr = real_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i; /*/(real_t)M;*/
  }
  
  j = real_stable_sum(ptr,N);

  k = (N/2.0)*(N-1.0); /* (N/(real_t)M) * ((N-(real_t)1.0)/(real_t)2.0); */

  TESTEQUALS(j,k,PF_REAL_T);

  return 0;
}
