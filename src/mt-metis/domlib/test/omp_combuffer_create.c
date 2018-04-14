#include "test.h"
#include <omp.h>

#define N 10

sint_t test(void) 
{
  int i;
  sint_t rv = 0;
  #pragma omp parallel shared(rv) num_threads(4)
  {
    const int myid = omp_get_thread_num();
    sint_combuffer_t * com = sint_combuffer_create(N);

    for (i=0;i<omp_get_num_threads();++i) {
      if (i != myid && com->buffers[omp_get_thread_num()][i].maxsize != N) {
        printf("[%d] buffer[%d]->maxsize = %zu and num threads = %d\n",
            myid,i,com->buffers[myid][i].maxsize,N);
        #pragma omp atomic
        ++rv;
      }
    }
    sint_combuffer_free(com);
  }
  return rv;
}
