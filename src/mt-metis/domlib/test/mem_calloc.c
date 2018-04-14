#include "test.h"

#define N 1000

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_calloc(N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],0,PF_SINT_T);
  }
  return 0;
}
