#include "test.h"

#define N 213523 

sint_t test(void)
{
  sint_t i,j;
  
  sint_t * perm = sint_alloc(N);
  sint_t * perm2 = sint_alloc(N);
  sint_t * mark = sint_calloc(N);

  sint_incset(perm,0,1,N);
  sint_incset(perm2,0,1,N);

  sint_shuffle(perm,N);
  sint_shuffle(perm2,N);

  for (i=0;i<N;++i) {
    j = perm[i];
    mark[j] = 1;
  }

  for (i=0;i<N;++i) {
    TESTEQUALS(mark[i],1,PF_SINT_T);
  }

  for (i=0;i<N;++i) {
    if (perm[i] != perm2[i]) {
      break;
    }
  }
  TESTLESSTHAN(i,N,PF_SINT_T);

  return 0;
}
