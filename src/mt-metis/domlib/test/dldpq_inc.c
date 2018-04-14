#include "test.h"

#define M 100
#define N 400
#define X -5
#define Y 150


sint_t test(void)
{
  sint_t i,j,k;
  
  sint_dldpq_t * q = sint_dldpq_create(M,N,X,Y);
  sint_t * keys = sint_init_alloc(M,Y-X);

  sint_dldpq_fill_min(q);

  for (i=X;i<Y;++i) {
    for (j=0;j<i;++j) {
      sint_dldpq_inc(i,q);
    }
  } 

  for (i=0;i<Y;++i) {
    j = sint_dldpq_peek(i,q); 
    TESTEQUALS(j,i+M,PF_SINT_T);
  }

  sint_dldpq_fill_min(q);

  for (i=X;i<Y;++i) {
    k = sint_rand(0,N-M);
    for (j=0;j<k;++j) {
      sint_dldpq_inc(i,q);
    }
    keys[i-X] = k+M;
  }

  for (i=X;i<Y;++i) {
    k = keys[i-X];
    j = sint_dldpq_peek(i,q);
    TESTEQUALS(j,k,PF_SINT_T);
  }

  sint_dldpq_fill_min(q);

  sint_set(keys,M,Y-X);
  for (i=M;i<N;++i) {
    k = sint_rand(X,Y);
    sint_dldpq_inc(k,q);
    ++keys[k-X];
  }

  for (i=X;i<Y;++i) {
    k = keys[i-X];
    j = sint_dldpq_peek(i,q);
    TESTEQUALS(j,k,PF_SINT_T);
  }

  sint_dldpq_free(q);

  return 0;
}
