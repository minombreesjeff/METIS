/**
 * @file dlomp.h
 * @brief OpenMP reduction functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-10-08
 */




#ifndef DL_OMP_H
#define DL_OMP_H




/******************************************************************************
* REDUCTION *******************************************************************
******************************************************************************/


#define __DL_MK_OMP_REDUCTION_HEADERS(prefix,type_t) \
  void prefix ## _omp_reduction_init(const size_t maxthreads); \
  void prefix ## _omp_reduction_clear(void); \
  type_t prefix ## _omp_maxreduce_value(type_t val); \
  size_t prefix ## _omp_maxreduce_index(type_t val); \
  type_t prefix ## _omp_minreduce_value(type_t val); \
  size_t prefix ## _omp_minreduce_index(type_t val); \
  type_t prefix ## _omp_sumreduce(type_t val); \
  void prefix ## _omp_sumareduce(type_t * val, size_t n);

#define __DL_MK_OMP_REDUCTION_FUNCS(prefix,type_t) \
  type_t * __ ## prefix ## _omp_redbuf = NULL; \
  type_t ** __ ## prefix ## _omp_dredbuf = NULL; \
  void prefix ## _omp_reduction_init(const size_t maxthreads) \
  { \
    _Pragma("omp barrier") \
    _Pragma("omp master") \
    { \
      __ ## prefix ## _omp_redbuf = prefix ## _alloc(maxthreads); \
      __ ## prefix ## _omp_dredbuf = \
        r_ ## prefix ## _alloc(size_uppow2(maxthreads)); \
    } \
    _Pragma("omp barrier") \
  } \
  void prefix ## _omp_reduction_clear(void) \
  { \
    _Pragma("omp barrier") \
    _Pragma("omp master") \
    { \
      dl_free(__ ## prefix ## _omp_redbuf); \
      dl_free(__ ## prefix ## _omp_dredbuf); \
      __ ## prefix ## _omp_redbuf = NULL; \
      __ ## prefix ## _omp_dredbuf = NULL; \
    } \
    _Pragma("omp barrier") \
  } \
  type_t prefix ## _omp_maxreduce_value(type_t val) \
  { \
    static type_t max; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    _Pragma("omp barrier") \
    __ ## prefix ## _omp_redbuf[myid] = val; \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      max = prefix ## _max_value(__ ## prefix ## _omp_redbuf, nthreads); \
    } \
    _Pragma("omp barrier") \
    return max; \
  } \
  size_t prefix ## _omp_maxreduce_index(type_t val) \
  { \
    static size_t max; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    _Pragma("omp barrier") \
    __ ## prefix ## _omp_redbuf[myid] = val; \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      max = prefix ## _max_index(__ ## prefix ## _omp_redbuf, nthreads); \
    } \
    _Pragma("omp barrier") \
    return max; \
  } \
  type_t prefix ## _omp_minreduce_value(type_t val) \
  { \
    static type_t min; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    _Pragma("omp barrier") \
    __ ## prefix ## _omp_redbuf[myid] = val; \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      min = prefix ## _min_value(__ ## prefix ## _omp_redbuf, nthreads); \
    } \
    _Pragma("omp barrier") \
    return min; \
  } \
  size_t prefix ## _omp_minreduce_index(type_t val) \
  { \
    static size_t min; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    _Pragma("omp barrier") \
    __ ## prefix ## _omp_redbuf[myid] = val; \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      min = prefix ## _min_index(__ ## prefix ## _omp_redbuf, nthreads); \
    } \
    _Pragma("omp barrier") \
    return min; \
  } \
  type_t prefix ## _omp_sumreduce(type_t val) \
  { \
    static type_t sum; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    _Pragma("omp barrier") \
    __ ## prefix ## _omp_redbuf[myid] = val; \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      sum = prefix ## _sum(__ ## prefix ## _omp_redbuf, nthreads); \
    } \
    _Pragma("omp barrier") \
    return sum; \
  } \
  static void __ ## prefix ## _omp_sumdreduce(type_t * const val, size_t n) \
  { \
    size_t i; \
    ssize_t d; \
    size_t nbr; \
    type_t * nbrval; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    const size_t up = size_uppow2(nthreads); \
    const size_t start = ((myid & 1U) == 0 ? 0 : n/2); \
    const size_t end = ((myid & 1U) == 0 ? n/2 : n); \
    const int sign = ((myid & 1U) == 0 ? 1 : -1);  \
    _Pragma("omp barrier") \
    __ ## prefix ## _omp_dredbuf[myid] = val; \
    if (myid + nthreads < up) { \
      __ ## prefix ## _omp_dredbuf[myid+nthreads] = prefix ## _calloc(n); \
    } \
    _Pragma("omp barrier") \
    nbr = myid; \
    for (d=1;d<(ssize_t)nthreads;d<<=1) { \
      nbr = (size_t)((nbr+(sign*d))%up); \
      nbrval = __ ## prefix ## _omp_dredbuf[nbr]; \
      if (nbr < nthreads) { \
        for (i=start;i<end;++i) { \
          val[i] += nbrval[i]; \
        } \
        prefix ## _copy(nbrval+start,val+start,end-start); \
      } else { \
        for (i=0;i<n;++i) { \
          val[i] += nbrval[i]; \
        } \
        prefix ## _copy(nbrval,val,n); \
      } \
      _Pragma("omp barrier") \
    } \
    if (myid + nthreads < up) { \
      dl_free(__ ## prefix ## _omp_dredbuf[myid+nthreads]); \
    } \
  } \
  void prefix ## _omp_sumareduce(type_t * const val, size_t n) \
  { \
    size_t i, j, t, ei; \
    const size_t chunk = 64 / sizeof(type_t); \
    type_t sum[chunk]; \
    const size_t myid = omp_get_thread_num(); \
    const size_t nthreads = omp_get_num_threads(); \
    if (nthreads > 128 && n / nthreads < chunk) { \
      __ ## prefix ## _omp_sumdreduce(val,n); \
    } else { \
      const size_t start = size ## _chunkstart(myid,nthreads,n); \
      const size_t end = start + size ## _chunksize(myid,nthreads,n); \
      __ ## prefix ## _omp_dredbuf[myid] = val; \
      _Pragma("omp barrier") \
      for (j=start;j<end;j+=chunk) { \
        ei = dl_min(chunk,end-j); \
        prefix ## _copy(sum,val+j,ei); \
        for (t=(myid+1)%nthreads;t!=myid;t=((t+1)%nthreads)) { \
          for (i=0;i<ei;++i) { \
            sum[i] += __ ## prefix ## _omp_dredbuf[t][j+i]; \
          } \
        } \
        prefix ## _copy(val+j,sum,ei); \
        for (t=(myid+1)%nthreads;t!=myid;t=((t+1)%nthreads)) { \
          prefix ## _copy(__ ## prefix ## _omp_dredbuf[t]+j,sum,ei); \
        } \
      } \
      _Pragma("omp barrier") \
    } \
  }


/* return the id of the thread with the minimum value */
#define dl_omp_prefixsum_reduce(myid,values,nvalues,partial,buffer,nthreads, \
    type) \
  do { \
    DL_ASSERT(nthreads>0,"Must have positive number of threads\n"); \
    _Pragma("omp barrier") \
    size_t _i; \
    size_t _mystart = ((nvalues/nthreads)*myid) + \
      dl_min((size_t)myid,(size_t)(nvalues%nthreads)); \
    size_t _mysize = (nvalues/nthreads) + \
      (((size_t)myid < (size_t)(nvalues%nthreads)) ? 1U : 0U); \
    for (_i=_mystart+1U;_i<_mystart+_mysize;++_i) { \
      values[_i] += values[_i-1U]; \
    } \
    partial[myid] = values[_i-1U]; \
    _Pragma("omp barrier") \
    for (_i=1U;_i<(size_t)nthreads;_i<<=1U) { \
      if (((size_t)myid)>=_i) { \
        buffer[myid] = partial[myid] + partial[myid-_i]; \
      } \
      _Pragma("omp barrier") \
      _Pragma("omp single") \
      { \
        memcpy(partial+1,buffer+1,sizeof(type)*(nthreads-1)); \
      } \
    } \
    for (_i=_mystart;_i<_mystart+_mysize;++_i) { \
      values[_i] += partial[myid] - values[_mystart+_mysize-1]; \
    } \
    _Pragma ("omp barrier") \
  } while(0)


/******************************************************************************
* RENAME MACROS ***************************************************************
******************************************************************************/
#define DL_MK_OMP_REDUCTION_HEADERS(prefix, type_t) \
  __DL_MK_OMP_REDUCTION_HEADERS(prefix, type_t)

#define DL_MK_OMP_REDUCTION_FUNCS(prefix, type_t) \
  __DL_MK_OMP_REDUCTION_FUNCS(prefix, type_t)

#endif
