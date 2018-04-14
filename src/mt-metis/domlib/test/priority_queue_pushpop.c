#include "test.h"

#define M 12000UL
#define N 16000UL




sint_t test(void)
{
  size_t l,i;
  real_t j;
  
  rs_priority_queue_t * q = rs_priority_queue_create(M,N);

  for (i=M;i<(sint_t)N;++i) {
    rs_maxpq_push(i/((real_t)N),i,q);
    for (l=0;l<q->size;++l) {
      TESTEQUALS(q->index[q->elements[l].val-q->min],l,PF_SIZE_T);
    }
  }

  TESTEQUALS(N-M,q->size,PF_SIZE_T);

  for (i=0;i<q->size;++i) {
    TESTEQUALS(q->index[((ssize_t)q->elements[i].val)-((ssize_t)q->min)],i,
        PF_SIZE_T);
  }

  j = q->elements[0].key;
  rs_maxpq_pop(q);
  for (i=1;i<N-M;++i) {
    TESTLESSTHANOREQUAL(q->elements[0].key,j,PF_REAL_T);
    j = q->elements[0].key;
    rs_maxpq_pop(q);
    for (l=0;l<q->size;++l) {
      TESTEQUALS(q->index[q->elements[l].val-q->min],l,PF_SIZE_T);
    }
  }

  TESTEQUALS(0UL,q->size,PF_SIZE_T);

  return 0;
}
