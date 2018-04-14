#include "test.h"

#define N 16000




sint_t test(void)
{
  sint_t i;
  real_t j,l;
  
  rs_priority_queue_t * q = rs_priority_queue_create(0,N);

  real_t * p = real_alloc(N);

  for (i=0;i<N;++i) {
    rs_maxpq_push((real_t)i,i,q);
    p[i] = (N%(i+10))/(real_t)N;
  }

  for (i=0;i<N;++i) {
    rs_maxpq_update(p[i],i,q);
  }

  l = rs_maxpq_max(q);
  while (q->size > 0) {
    i = rs_maxpq_peek(q);
    j = rs_maxpq_max(q);
    TESTEQUALS(p[i],j,PF_REAL_T);
    rs_maxpq_pop(q);
    /* make sure they priority queue still works */
    TESTLESSTHANOREQUAL(j,l,PF_REAL_T);
  }

  dl_free(p);
  rs_priority_queue_free(q);

  return 0;
}
