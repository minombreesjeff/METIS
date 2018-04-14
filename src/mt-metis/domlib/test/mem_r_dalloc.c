#include "test.h"

#define N 5

sint_t test(void) 
{
  sint_t i,j;
  sint_t ns[] = {1,2,3,4,5};
  sint_t ** ptr = r_sint_dalloc(ns,sizeof(sint_t),N);

  for (i=0;i<N;++i) {
    TESTTRUE(ptr[i]!=NULL);
    for (j=0;j<ns[i];++j) {
      ptr[i][j] = 0;
      TESTEQUALS(ptr[i][j],0,PF_SINT_T);
    }
    dl_free(ptr[i]);
  }
  
  dl_free(ptr);
  return 0;
}
