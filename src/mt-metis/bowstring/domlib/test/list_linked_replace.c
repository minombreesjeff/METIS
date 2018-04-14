#include "test.h"

#define N (10000)

sint_t test(void)
{
  sint_t i,j;
  sint_list_t * list = sint_list_create();

  for (i=0;i<N;++i) {
    sint_list_add(i,list);
  }

  for (i=0;i<N;++i) {
    j = sint_list_replace(N-i,i,list);
    TESTEQUALS(j,i,PF_SINT_T);
  }

  for (i=0;i<N;++i) {
    TESTEQUALS(sint_list_get(i,list),N-i,PF_SINT_T);
  }

  sint_list_free(list);

  return 0;
}
