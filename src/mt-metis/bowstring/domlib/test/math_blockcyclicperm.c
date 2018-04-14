#include "test.h"

#define N (10)
#define M (12345678)
#define O (1382908)

sint_t test(void) 
{
  sint_t i;
  sint_t * ptr = sint_alloc(N);

  sint_t ans1[] = {0,1,6,7,2,3,8,9,4,5};
  sint_t ans2[] = {0,1,2,6,7,8,3,4,5,9};
  sint_t ans3[] = {0,1,2,3,4,5,6,7,8,9};
  sint_t ans4[] = {0,1,2,3,4,5,6,7,8,9};
  sint_t ans5[] = {0,3,6,9,1,4,7,2,5,8};

  sint_blockcyclicperm(ptr,3,2,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans1[i],PF_SINT_T);
  }

  sint_blockcyclicperm(ptr,2,3,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans2[i],PF_SINT_T);
  }

  sint_blockcyclicperm(ptr,3,5,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans3[i],PF_SINT_T);
  }

  sint_blockcyclicperm(ptr,4,4,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans4[i],PF_SINT_T);
  }

  sint_blockcyclicperm(ptr,3,1,N);

  for (i=0;i<N;++i) {
    TESTEQUALS(ptr[i],ans5[i],PF_SINT_T);
  }

  dl_free(ptr);

  ptr = sint_init_alloc(-1,M);
  sint_blockcyclicperm(ptr,37,219,M);

  for (i=0;i<M;++i) {
    TESTGREATERTHANOREQUAL(ptr[i],0,PF_SINT_T);
    TESTLESSTHAN(ptr[i],M,PF_SINT_T);
  }

  dl_free(ptr);

  ptr = sint_init_alloc(-1,O);
  sint_blockcyclicperm(ptr,4,4096,O);

  for (i=0;i<O;++i) {
    TESTGREATERTHANOREQUAL(ptr[i],0,PF_SINT_T);
    TESTLESSTHAN(ptr[i],O,PF_SINT_T);
  }

  dl_free(ptr);



  return 0;
}
