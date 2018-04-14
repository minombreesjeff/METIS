#include "test.h"

#define M 500UL
#define N 5000UL




sint_t test(void) 
{
  sint_t i,j;

  sint_buffer_t * buffer = sint_buffer_create(M);

  TESTEQUALS(0UL,buffer->size,PF_SIZE_T);

  for (i=0;i<N;++i) {
    sint_buffer_add(i,buffer);
  }

  TESTEQUALS(N,buffer->size,PF_SIZE_T);

  for (i=0;i<N;++i) {
    TESTEQUALS(i,buffer->elements[i],PF_SINT_T);
  }

  sint_buffer_free(buffer);

  return 0;
}
