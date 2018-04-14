#include "test.h"

#define K 120500UL
#define L 121024UL
#define M 120000UL
#define N 160000UL




sint_t test(void)
{
  size_t i;
  
  rs_priority_queue_t * q = rs_priority_queue_create(M,N);

  for (i=M;i<N;++i) {
    rs_maxpq_push(i/((real_t)N),i,q);
  }

  TESTEQUALS(N-M,q->size,PF_SIZE_T);

  rs_maxpq_clear(q);

  for (i=0;i<q->maxsize;++i) {
    TESTEQUALS(q->index[i],(size_t)-1,PF_SIZE_T);
  }

  rs_priority_queue_free(q);

  return 0;
}
