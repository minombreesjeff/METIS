#include "test.h"

#define K 120500UL
#define L 121024UL
#define M 120000UL
#define N 160000UL




sint_t test(void)
{
  size_t i;
  real_t j;
  
  rs_priority_queue_t * q = rs_priority_queue_create(M,N);

  for (i=M;i<N;++i) {
    rs_maxpq_push(i/((real_t)N),i,q);
  }

  TESTEQUALS(N-M,q->size,PF_SIZE_T);

  for (i=K;i<L;++i) {
    rs_maxpq_delete(i,q);
  }

  TESTEQUALS(N-M-(L-K),q->size,PF_SIZE_T);

  j = q->elements[0].key;
  rs_maxpq_pop(q);
  while (q->size > 0) {
    TESTLESSTHANOREQUAL(q->elements[0].key,j,PF_REAL_T);
    j = q->elements[0].key;
    rs_maxpq_pop(q); 
  }

  for (i=0;i<q->maxsize;++i) {
    TESTEQUALS(q->index[i],(size_t)-1,PF_SIZE_T);
  }

  rs_priority_queue_free(q);

  return 0;
}
