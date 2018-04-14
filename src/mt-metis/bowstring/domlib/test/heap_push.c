#include "test.h"

#define N 160000UL
#define M 137UL

sint_t test(void)
{
  sint_t i;

  sint_heap_t * maxh = sint_heap_create(M);
  
  for (i=0;i<(sint_t)N;++i) {
    sint_heap_push(i,maxh);
  }

  TESTEQUALS((sint_t)(N-1),maxh->elements[0],PF_SINT_T);

  for (i=1;i<(sint_t)N;++i) {
    TESTGREATERTHANOREQUAL(maxh->elements[(i-1)>>1],maxh->elements[i],
        PF_SINT_T);
  }

  TESTEQUALS(maxh->size,N,PF_SIZE_T);

  sint_heap_free(maxh);

  return 0;
}
