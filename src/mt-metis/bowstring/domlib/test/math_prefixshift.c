#include "test.h"

#define N 1000

sint_t test(void) 
{
  sint_t i, end;
  sint_t * ptr = sint_alloc(N);
  sint_t * shift = sint_alloc(N+1);

  for (i=0;i<N;++i) {
    ptr[i] = sint_rand(0,100);
  }

  shift[0] = 0;
  sint_copy(shift+1,ptr,N);

  end = sint_prefixshift(ptr,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(shift[i],ptr[i],PF_SINT_T);
  }
  TESTEQUALS(end,shift[N],PF_SINT_T);
  return 0;
}
