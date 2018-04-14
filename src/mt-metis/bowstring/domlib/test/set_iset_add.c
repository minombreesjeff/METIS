#include "test.h"

#define M 500UL
#define N 1000UL

sint_t test(void) 
{
  sint_t i,j;

  sint_iset_t * set = sint_iset_create(M,N);

  TESTEQUALS(0UL,set->size,PF_SIZE_T);

  for (i=M;i<(sint_t)N;++i) {
    if (i % 2 == 0) {
      sint_iset_add(i,set);
    }
  }

  TESTEQUALS((N-M)/2UL,set->size,PF_SIZE_T);

  j=0;
  for (i=M;i<(sint_t)N;++i) {
    if (i % 2 == 0) {
      TESTEQUALS(i,set->ind[j],PF_SINT_T);
      ++j;
    } else {
      TESTEQUALS(-1,set->ptr[i-M],PF_SINT_T);
    }
  }

  sint_iset_free(set);

  return 0;
}
