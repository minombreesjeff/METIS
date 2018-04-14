#include "test.h"

#define N 160000

/* I should re-write this */

sint_t test(void)
{
  unsigned int seed = 1;
  sint_t i;
  
  sint_t * unsorted = sint_alloc(N);
  sint_t * unsortedv = sint_alloc(N);
  sint_t * sorted = sint_alloc(N);

  sint_incset(unsorted,0,1,N);
  sint_incset(unsortedv,0,1,N);

  sint_shuffle_r(unsorted,N,&seed);
  seed = 1;
  sint_shuffle_r(unsortedv,N,&seed);

  sint_countingsort_v(unsorted,unsortedv,sorted,0,N-1,N);

  sint_incset(unsorted,0,1,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(unsorted[i],sorted[i],PF_SINT_T);
  }

  dl_free(unsorted);
  dl_free(unsortedv);
  dl_free(sorted);

  return 0;
}
