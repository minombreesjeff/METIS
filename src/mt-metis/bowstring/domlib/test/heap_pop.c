#include "test.h"

#define N 160000
#define M 137

sint_t test(void)
{
  sint_t i,j;

  sint_heap_t * maxh = sint_heap_create(N);
  

  for (i=1;(i<<1)<N;i<<=1) {
    j = (i<<1)-1;
    sint_incset(maxh->elements+i-1,N-j,1,i);
  }
  sint_incset(maxh->elements+i-1,0,1,N-i+1);


  maxh->size = N;

  for (i=0;i<N;++i) {
    j = sint_heap_pop(maxh);
    TESTEQUALS(N-i-1,j,PF_SINT_T);
  }

  TESTEQUALS(maxh->size,0UL,PF_SIZE_T);

  sint_heap_free(maxh);

  return 0;
}
