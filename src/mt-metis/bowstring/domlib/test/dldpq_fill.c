#include "test.h"

#define M 100
#define N 400
#define X -5
#define Y 150


sint_t test(void)
{
  sint_t i,j;
  
  sint_dldpq_t * q = sint_dldpq_create(M,N,X,Y);

  sint_dldpq_fill_min(q);

  for (i=X;i<Y;++i) {
    j = sint_dldpq_peek(i,q);
    TESTEQUALS(j,M,PF_SINT_T);
  } 

  sint_dldpq_fill_max(q);

  for (i=X;i<Y;++i) {
    j = sint_dldpq_peek(i,q);
    TESTEQUALS(j,N-1,PF_SINT_T);
  } 

  sint_dldpq_fill_min_perm(q);

  for (i=X;i<Y;++i) {
    j = sint_dldpq_peek(i,q);
    TESTEQUALS(j,M,PF_SINT_T);
  } 

  sint_dldpq_fill_max_perm(q);

  for (i=X;i<Y;++i) {
    j = sint_dldpq_peek(i,q);
    TESTEQUALS(j,N-1,PF_SINT_T);
  } 

  sint_dldpq_free(q);

  return 0;
}
