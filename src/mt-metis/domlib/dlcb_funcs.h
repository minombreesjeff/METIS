/**
 * @file dlcb_funcs.h
 * @brief Functions for communication buffers
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-10-05
 */




/* prefixing ugliness */
#define DLCB_PRE2(prefix,suffix) prefix ## _ ## suffix
#define DLCB_PRE1(prefix,suffix) DLCB_PRE2(prefix,suffix)
#define DLCB_PUB(name) DLCB_PRE1(DLCB_PREFIX,name)
#define DLCB_PRI(name) DLCB_PRE1(_,DLCB_PRE1(DLCB_PREFIX,name))
#define DLCB_RPRI(name) DLCB_PRE1(r__,DLCB_PRE1(DLCB_PREFIX,name))




/******************************************************************************
* INCLUDES ********************************************************************
******************************************************************************/


#ifndef DLCB_BUFFER_TYPE_T

#define DLBUFFER_PREFIX DLCB_PRI(cb)
#define DLBUFFER_TYPE_T DLCB_TYPE_T
#include "dlbuffer_funcs.h"
#undef DLBUFFER_PREFIX
#undef DLBUFFER_TYPE_T

#define DLCB_BUFFER_FUNC(func) DLCB_PRE1(DLCB_PRI(cb_buffer),func)

#else

#define DLCB_BUFFER_FUNC(func) DLCB_PRE1(DLCB_BUFFER_PREFIX,func)

#endif


#define DLMEM_PREFIX DLCB_PUB(combuffer)
#define DLMEM_TYPE_T DLCB_PUB(combuffer_t)
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


#ifndef DLCB_VISIBILITY
  #define DLCB_DEFVIS
  #define DLCB_VISIBILITY
#endif




DLCB_VISIBILITY DLCB_PUB(combuffer_t) * DLCB_PRI(cbc_com);
DLCB_VISIBILITY DLCB_PUB(combuffer_t) * DLCB_PUB(combuffer_create)(
    const size_t n)
{
  size_t i;
  const size_t myid = omp_get_thread_num();
  #pragma omp barrier
  #pragma omp master
  {
    DLCB_PRI(cbc_com) = DLCB_PUB(combuffer_alloc)(1);
    DLCB_PRI(cbc_com)->nthreads = omp_get_num_threads();
    #ifndef DLCB_BUFFER_TYPE_T
    DLCB_PRI(cbc_com)->buffers = (DLCB_PRI(cb_buffer_t)**)
    #else
    DLCB_PRI(cbc_com)->buffers = (DLCB_BUFFER_TYPE_T**)
    #endif
        malloc(sizeof(void*)*DLCB_PRI(cbc_com)->nthreads);
  }
  #pragma omp barrier
  DLCB_PRI(cbc_com)->buffers[myid] =
      DLCB_BUFFER_FUNC(alloc)(DLCB_PRI(cbc_com)->nthreads);
  for (i=0;i<DLCB_PRI(cbc_com)->nthreads;++i) {
    if (i != myid) { 
      DLCB_BUFFER_FUNC(init)(n,DLCB_PRI(cbc_com)->buffers[myid]+i);
    } else {
      DLCB_PRI(cbc_com)->buffers[myid][i].maxsize = 0;
      DLCB_PRI(cbc_com)->buffers[myid][i].defsize = 0;
      DLCB_PRI(cbc_com)->buffers[myid][i].size = 0;
      DLCB_PRI(cbc_com)->buffers[myid][i].elements = NULL;
    }
  }
  return DLCB_PRI(cbc_com);
}


DLCB_VISIBILITY void DLCB_PUB(combuffer_add)(const size_t dst, 
    const DLCB_TYPE_T val, DLCB_PUB(combuffer_t) * const com)
{
  DLCB_BUFFER_FUNC(add)(val,com->buffers[omp_get_thread_num()]+dst);
}


DLCB_VISIBILITY void DLCB_PUB(combuffer_send)(DLCB_PUB(combuffer_t) * com)
{
  size_t i,j;
  DLCB_BUFFER_FUNC(t) buf;
  #pragma omp barrier
  #pragma omp for schedule(static)
  for (i=0;i<com->nthreads-1;++i) {
    for (j=i+1;j<com->nthreads;++j) {
      buf = com->buffers[j][i];
      com->buffers[j][i] = com->buffers[i][j];
      com->buffers[i][j] = buf;
    }
  }
  #pragma omp barrier
}


DLCB_VISIBILITY void DLCB_PUB(combuffer_clear)(DLCB_PUB(combuffer_t) * com)
{
  size_t i;
  const size_t myid = omp_get_thread_num();
  for (i=0;i<com->nthreads;++i) {
    DLCB_BUFFER_FUNC(clear)(com->buffers[myid]+i);
  }
}


DLCB_VISIBILITY void DLCB_PUB(combuffer_free)(DLCB_PUB(combuffer_t) * com)
{
  size_t i;
  const size_t myid = omp_get_thread_num();
  for (i=0;i<com->nthreads;++i) {
    if (i!=myid) {
      dl_free(com->buffers[myid][i].elements);
    }
  }
  dl_free(com->buffers[myid]);
  #pragma omp barrier
  #pragma omp master
  {
    dl_free(com->buffers);
    dl_free(com);
  }
}


#ifdef DLCB_BUFFER_TYPE_T

DLCB_VISIBILITY DLCB_BUFFER_TYPE_T * DLCB_PUB(combuffer_get)(const size_t idx, 
    const DLCB_PUB(combuffer_t) * const com)
{
  return com->buffers[omp_get_thread_num()]+idx;
}

#endif



#ifdef DLCB_DEFVIS
  #undef DLCB_DEFVIS
  #undef DLCB_VISIBILITY
#endif



#undef DLCB_PRE
#undef DLCB_PRE
#undef DLCB_PUB
#undef DLCB_PRI
#undef DLCB_RPRI
#undef DLCB_BUFFER_FUNC



