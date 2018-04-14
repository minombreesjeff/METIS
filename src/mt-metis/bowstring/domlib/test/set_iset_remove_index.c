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

  /* random remove */
  for (i=M;i<N;++i) {
    __LSTADD(i,set);
  }
  while (set->size > 0) {
    i = sint_rand(0,set->size);
    j = sint_iset_remove_index(i,set);
    TESTGREATERTHANOREQUAL(j,M,PF_SINT_T);
    TESTLESSTHAN(j,N,PF_SINT_T);
  }

  TESTEQUALS(0UL,set->size,PF_SIZE_T);

  /* remove from end */
  for (i=M;i<N;++i) {
    __LSTADD(i,set);
  }
  while (set->size > 0) {
    j = sint_iset_remove_index(set->size-1,set);
    TESTGREATERTHANOREQUAL(j,M,PF_SINT_T);
    TESTLESSTHAN(j,N,PF_SINT_T);
  }

  TESTEQUALS(0UL,set->size,PF_SIZE_T);

  sint_iset_free(set);

  return 0;
}
