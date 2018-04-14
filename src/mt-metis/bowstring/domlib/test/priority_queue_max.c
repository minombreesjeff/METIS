#include "test.h"

#define N 1600UL




sint_t test(void)
{
  size_t i;
  real_t j;
  
  rs_priority_queue_t * q = rs_priority_queue_create(0,N);

  real_t * p = real_alloc(N);

  for (i=0;i<N;++i) {
    p[i] = (N % (i+10))/(real_t)N;
    rs_maxpq_push(p[i],i,q);
  }

  while (q->size > 0) {
    j = rs_maxpq_max(q);
    i = rs_maxpq_pop(q);
    TESTEQUALS(p[i],j,PF_REAL_T);
  }

  dl_free(p);
  rs_priority_queue_free(q);

  return 0;
}
