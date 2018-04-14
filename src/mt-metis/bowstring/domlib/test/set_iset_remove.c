#include "test.h"

#define M 500
#define N 1000

#define __LSTADD(i,lst) \
  do { \
    (lst)->ind[(lst)->size] = i; \
    (lst)->ptr[i-(lst)->min] = (lst)->size++; \
  } while (0) 




sint_t test(void) 
{
  sint_t i,j;

  sint_iset_t * set = sint_iset_create(M,N);

  TESTEQUALS(0UL,set->size,PF_SIZE_T);

  for (i=M;i<N;++i) {
    __LSTADD(i,set);
  }

  for (i=M;i<N;++i) {
    if (i % 2 == 1) {
      sint_iset_remove(i,set);
    }
  }

  j=0;
  for (i=M;i<N;++i) {
    if (i % 2 == 0) {
      TESTEQUALS(j,set->ptr[set->ind[j]-set->min],PF_SINT_T);
      ++j;
    } else {
      TESTEQUALS(-1,set->ptr[i-M],PF_SINT_T);
    }
  }

  sint_iset_free(set);

  return 0;
}
