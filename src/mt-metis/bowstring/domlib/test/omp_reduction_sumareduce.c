#include "test.h"
#include <omp.h>

#define N 1000000

sint_t test(void) 
{
  sint_t ans;
  sint_t rv = 0;
  sint_omp_reduction_init(64);
  ans = 0;
  #pragma omp parallel shared(rv,ans) num_threads(1)
  {
    sint_t i;
    sint_t * mybuf = sint_init_alloc(omp_get_thread_num(),N);
    sint_omp_sumareduce(mybuf,N);
    for (i=0;i<N;++i) {
      OMPTESTEQUALS(ans,mybuf[i],PF_SINT_T,rv);
    }
    dl_free(mybuf);
  }
  if (rv != 0) {
    return rv;
  }
  ans=3;
  #pragma omp parallel shared(rv,ans) num_threads(3)
  {
    sint_t i;
    sint_t * mybuf = sint_init_alloc(omp_get_thread_num(),N);
    sint_omp_sumareduce(mybuf,N);
    for (i=0;i<N;++i) {
      OMPTESTEQUALS(ans,mybuf[i],PF_SINT_T,rv);
    }
    dl_free(mybuf);
  }
  if (rv != 0) {
    return rv;
  }
  ans=6;
  #pragma omp parallel shared(rv,ans) num_threads(4)
  {
    sint_t i;
    sint_t * mybuf = sint_init_alloc(omp_get_thread_num(),N);
    sint_omp_sumareduce(mybuf,N);
    for (i=0;i<N;++i) {
      OMPTESTEQUALS(ans,mybuf[i],PF_SINT_T,rv);
    }
    dl_free(mybuf);
  }
  if (rv != 0) {
    return rv;
  }
  ans=28;
  #pragma omp parallel shared(rv,ans) num_threads(8)
  {
    sint_t i;
    sint_t * mybuf = sint_init_alloc(omp_get_thread_num(),N);
    sint_omp_sumareduce(mybuf,N);
    for (i=0;i<N;++i) {
      OMPTESTEQUALS(ans,mybuf[i],PF_SINT_T,rv);
    }
    dl_free(mybuf);
  }
  if (rv != 0) {
    return rv;
  }
  ans=253;
  #pragma omp parallel shared(rv,ans) num_threads(23)
  {
    sint_t i;
    sint_t * mybuf = sint_init_alloc(omp_get_thread_num(),N);
    sint_omp_sumareduce(mybuf,N);
    for (i=0;i<N;++i) {
      OMPTESTEQUALS(ans,mybuf[i],PF_SINT_T,rv);
    }
    dl_free(mybuf);
  }
  if (rv != 0) {
    return rv;
  }
  ans=2016;
  #pragma omp parallel shared(rv,ans) num_threads(64)
  {
    sint_t i;
    sint_t * mybuf = sint_init_alloc(omp_get_thread_num(),N);
    sint_omp_sumareduce(mybuf,N);
    for (i=0;i<N;++i) {
      OMPTESTEQUALS(ans,mybuf[i],PF_SINT_T,rv);
    }
    dl_free(mybuf);
  }
  sint_omp_reduction_clear();
  return rv;
}
