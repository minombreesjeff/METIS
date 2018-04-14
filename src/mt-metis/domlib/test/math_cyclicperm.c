#include "test.h"

#define N (10)

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  sint_t ans1[] = {0,7,1,8,2,9,3,4,5,6};
  sint_t ans2[] = {0,3,6,9,1,4,7,2,5,8};

  sint_cyclicperm(ptr,7,10);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans1[i],PF_SINT_T);
  }

  sint_cyclicperm(ptr,3,10);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans2[i],PF_SINT_T);
  }

  dl_free(ptr);

  return 0;
}
