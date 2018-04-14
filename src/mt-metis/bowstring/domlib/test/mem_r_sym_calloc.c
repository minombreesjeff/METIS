#include "test.h"

#define N 1000
#define M 100

sint_t test(void) 
{
  sint_t i,j;
  sint_t ** ptr = r_sint_sym_calloc(M,N);

  for (i=0;i<N;++i) {
    for (j=0;j<M;++j) {
      TESTEQUALS(ptr[i][j],0,PF_SINT_T);
    }
    dl_free(ptr[i]);
  }

  dl_free(ptr);

  return 0;
}
