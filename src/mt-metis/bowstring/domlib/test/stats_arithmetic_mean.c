#include "test.h"

#define N 1000000

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  double emean = 0.0;
  double amean; 

  for (i=0;i<N;++i) {
    ptr[i] = sint_rand(0,100);
  }

  for (i=0;i<N;++i) {
    emean += ptr[i];
  }
  emean /= N;

  amean = sint_arithmetic_mean(ptr,N);

  TESTEQUALS(emean,amean,"%lf");

  return 0;
}
