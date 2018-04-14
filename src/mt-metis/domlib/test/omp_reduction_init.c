#include "test.h"
#include <omp.h>

#define N 10

extern sint_t * __sint_omp_redbuf;
extern sint_t ** __sint_omp_dredbuf;

sint_t test(void) 
{
  sint_t rv = 0;
  sint_omp_reduction_init(omp_get_num_threads());
  TESTTRUE((void*)__sint_omp_redbuf != NULL);
  TESTTRUE((void*)__sint_omp_dredbuf != NULL);
  sint_omp_reduction_clear();
  TESTEQUALS((void*)__sint_omp_redbuf,NULL,"%p");
  TESTEQUALS((void*)__sint_omp_dredbuf,NULL,"%p");
  #pragma omp parallel shared(rv) num_threads(4)
  {
    sint_omp_reduction_init(omp_get_num_threads());
    OMPTESTTRUE((void*)__sint_omp_redbuf != NULL,rv);
    OMPTESTTRUE((void*)__sint_omp_dredbuf != NULL,rv);
    sint_omp_reduction_clear();
  }
  TESTEQUALS((void*)__sint_omp_redbuf,NULL,"%p");
  TESTEQUALS((void*)__sint_omp_dredbuf,NULL,"%p");
  return rv;
}
