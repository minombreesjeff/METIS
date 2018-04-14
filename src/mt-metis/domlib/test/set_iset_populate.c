#include "test.h"

#define L 0UL
#define M 500UL
#define N 1000UL

sint_t test(void) 
{
  sint_t i;

  sint_iset_t * set = sint_iset_create(M,N);
  sint_iset_t * set2 = sint_iset_create(M,N);

  TESTEQUALS(0UL,set->size,PF_SIZE_T);

  sint_iset_populate(set);

  TESTEQUALS((N-M),set->size,PF_SIZE_T);

  for (i=M;i<(sint_t)N;++i) {
    sint_iset_add(i,set2);
  }

  for (i=M;i<(sint_t)N;++i) {
    TESTEQUALS((sint_t)(i-M),sint_iset_indexof(i,set),PF_SINT_T);
    TESTEQUALS(i,sint_iset_get(i-M,set),PF_SINT_T);
    TESTEQUALS(sint_iset_get(i-M,set2),sint_iset_get(i-M,set),PF_SINT_T);
    TESTTRUE(sint_iset_contains(i,set));
  }


  sint_iset_free(set);
  sint_iset_free(set2);

  set = sint_iset_create(L,N);
  set2 = sint_iset_create(L,N);

  sint_iset_populate(set);

  for (i=L;i<(sint_t)N;++i) {
    sint_iset_add(i,set2);
  }

  for (i=L;i<(sint_t)N;++i) {
    TESTEQUALS((sint_t)(i-L),sint_iset_indexof(i,set),PF_SINT_T);
    TESTEQUALS(i,sint_iset_get(i-L,set),PF_SINT_T);
    TESTEQUALS(sint_iset_get(i-L,set2),sint_iset_get(i-L,set),PF_SINT_T);
    TESTTRUE(sint_iset_contains(i,set));
  }

  sint_iset_free(set);
  sint_iset_free(set2);

  return 0;
}
