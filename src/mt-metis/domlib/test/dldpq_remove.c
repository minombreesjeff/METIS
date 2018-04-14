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
    for (j=X;j<i;++j) {
      sint_dldpq_inc(i,q);
    }
  } 

  k = Y - 1;
  for (i=X;i<Y;++i) {
    j = sint_dldpq_remove_max(q); 
    TESTEQUALS(j,k,PF_SINT_T);
    --k;
  }

  sint_dldpq_fill_min(q);

  for (i=X;i<Y;++i) {
    if (i%2 == 0) {
      for (j=X;j<i;++j) {
        sint_dldpq_inc(i,q);
      }
    }
  } 

  k = Y - 2;
  for (i=0;i<((Y-X)/2);++i) {
    j = sint_dldpq_remove_max(q); 
    TESTEQUALS(j,k,PF_SINT_T);
    k -= 2;
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
    j = sint_dldpq_remove(i,q);
    TESTEQUALS(j,k,PF_SINT_T);
  }

  sint_dldpq_fill_min(q);

  for (i=X;i<Y;++i) {
    k = sint_rand(0,N-M);
    for (j=0;j<k;++j) {
      sint_dldpq_inc(i,q);
    }
  }

  k = N;
  for (i=X;i<Y;++i) {
    j = sint_dldpq_peek_max(q);
    j = sint_dldpq_remove(j,q);
    TESTLESSTHANOREQUAL(j,k,PF_SINT_T);
    k = j;
  }

  sint_dldpq_free(q);

  return 0;
}
