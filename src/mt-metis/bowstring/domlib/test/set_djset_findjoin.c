#include "test.h"

#define M 500UL
#define N 1000UL

sint_t test(void) 
{
  sint_t i;

  sint_djset_t * set = sint_djset_create(M,N);

  TESTEQUALS(N-M,set->nsets,PF_SIZE_T);

  for (i=(sint_t)M;i<(sint_t)(N-1);++i) {
    sint_djset_join(i,N-1,set);
  }

  TESTEQUALS(sint_djset_find(M,set),sint_djset_find(N-5,set),PF_SINT_T);
  TESTEQUALS(1UL,set->nsets,PF_SIZE_T);

  sint_djset_free(set);

  set = sint_djset_create(0,10); 

  sint_djset_join(0,1,set);
  sint_djset_join(2,3,set);
  sint_djset_join(7,9,set);
  sint_djset_join(5,8,set);

  TESTEQUALS(6UL,set->nsets,PF_SIZE_T); 

  sint_djset_join(0,4,set);
  sint_djset_join(4,6,set);

  TESTEQUALS(4UL,set->nsets,PF_SIZE_T); 

  sint_djset_join(9,1,set);
  sint_djset_join(7,2,set);

  TESTEQUALS(2UL,set->nsets,PF_SIZE_T); 

  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(1,set),PF_SINT_T);
  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(2,set),PF_SINT_T);
  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(3,set),PF_SINT_T);
  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(4,set),PF_SINT_T);
  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(6,set),PF_SINT_T);
  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(7,set),PF_SINT_T);
  TESTEQUALS(sint_djset_find(0,set),sint_djset_find(9,set),PF_SINT_T);

  TESTEQUALS(sint_djset_find(5,set),sint_djset_find(8,set),PF_SINT_T);

  TESTTRUE(sint_djset_find(0,set) != sint_djset_find(8,set));

  sint_djset_free(set);

  return 0;
}
