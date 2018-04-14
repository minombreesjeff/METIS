#include "test.h"

#define M 500UL
#define N 5000UL




sint_t test(void) 
{
  sint_t i,j;

  sint_buffer_t * buffer = sint_buffer_create(M);

  TESTEQUALS(0UL,buffer->size,PF_SIZE_T);

  buffer->elements = sint_realloc(buffer->elements,N);
  buffer->maxsize = N;
  buffer->size = N;
  sint_incset(buffer->elements,N);

  sint_buffer_clear(buffer);

  TESTEQUALS(0,buffer->size,PF_SIZE_T);
  TESTEQUALS(N,buffer->maxsize,PF_SIZE_T);

  sint_incset(buffer->elements,M);
  buffer->size = M;

  sint_buffer_reset(buffer);

  TESTEQUALS(0,buffer->size,PF_SIZE_T);
  TESTEQUALS(M,buffer->maxsize,PF_SIZE_T);

  sint_buffer_free(buffer);

  return 0;
}
