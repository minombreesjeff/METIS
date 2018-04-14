#include "test.h"
#include <omp.h>

#define N 10UL
#define M 200UL
#define NTHREADS 64UL






sint_t test(void) {
  sint_t rv = 0;
  #pragma omp parallel shared(rv) num_threads(NTHREADS)
  {
    int i,j;
    const int myid = omp_get_thread_num();
    sint_combuffer_t * com = sint_combuffer_create(N);

    for (i=0;i<omp_get_num_threads();++i) {
      if (i != myid) {
        for (j=0;j<(int)M;++j) {
          sint_combuffer_add(i,myid,com);
        }
        OMPTESTEQUALS(M,com->buffers[myid][i].size,PF_SIZE_T,rv);
      }
    }
    for (i=0;i<omp_get_num_threads();++i) {
      if (i != myid) {
        for (j=0;j<(int)M;++j) {
          OMPTESTEQUALS(myid,com->buffers[myid][i].elements[j],"%d",rv);
        }
      }
    }
    sint_combuffer_free(com);
  }
  return rv;
}
