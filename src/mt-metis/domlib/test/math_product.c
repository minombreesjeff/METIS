#include "test.h"

#define N 100

sint_t test(void) 
{
  sint_t i,j,k;
  sint_t * ptr = sint_alloc(N);

  for (i=0;i<N;++i) {
    ptr[i] = i%5;
  }
  
  j = sint_product(ptr,N);

  k = 1;
  for (i=0;i<N;++i) {
    k *= ptr[i]; 
  }

  TESTEQUALS(j,k,PF_SINT_T);
  return 0;
}
