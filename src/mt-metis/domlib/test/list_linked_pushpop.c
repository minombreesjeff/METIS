#include "test.h"

#define N (10000)

sint_t test(void)
{
  sint_t i,j;
  sint_list_t * list = sint_list_create();

  for (i=0;i<N;++i) {
    sint_list_push(i,list);
    TESTEQUALS(list->size,i+1UL,"%zu");
  }

  for (i=N-1;i>=0;--i) {
    j = sint_list_pop(list);
    TESTEQUALS(j,i,PF_SINT_T);
    TESTEQUALS(list->size,(size_t)i,"%zu");
  }

  sint_list_free(list);

  return 0;
}
