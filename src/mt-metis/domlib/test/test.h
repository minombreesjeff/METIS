#ifndef TEST_H
#define TEST_H

#include <domlib.h>
#include <omp.h>


typedef float real_t;
#define PF_REAL_T "%f"
typedef int sint_t;
#define PF_SINT_T "%d"
typedef unsigned int uint_t;
#define PF_UINT_T "%u"


/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX real
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLPQ_PREFIX rs
#define DLPQ_KEY_T real_t
#define DLPQ_VAL_T sint_t
#include "dlpq_headers.h"
#undef DLPQ_PREFIX
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T




/* sint_t */
#define DLMEM_PREFIX sint 
#define DLMEM_TYPE_T sint_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX sint
#define DLMATH_TYPE_T sint_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLBUFFER_PREFIX sint
#define DLBUFFER_TYPE_T sint_t
#include "dlbuffer_headers.h"
#undef DLBUFFER_PREFIX
#undef DLBUFFER_TYPE_T


#define DLHT_PREFIX ss
#define DLHT_KEY_T sint_t
#define DLHT_VAL_T sint_t
#include "dlht_headers.h"
#undef DLHT_PREFIX
#undef DLHT_KEY_T
#undef DLHT_VAL_T


#define DLISET_PREFIX sint
#define DLISET_TYPE_T sint_t
#include "dliset_headers.h"
#undef DLISET_PREFIX
#undef DLISET_TYPE_T


#define DLDJSET_PREFIX sint
#define DLDJSET_TYPE_T sint_t
#include "dldjset_headers.h"
#undef DLDJSET_PREFIX
#undef DLDJSET_TYPE_T


#define DLCB_PREFIX sint
#define DLCB_TYPE_T sint_t
#include "dlcb_headers.h"
#undef DLCB_PREFIX
#undef DLCB_TYPE_T


#define DLDPQ_PREFIX sint
#define DLDPQ_KEY_T sint_t
#define DLDPQ_VAL_T sint_t
#include "dldpq_headers.h"
#undef DLDPQ_PREFIX
#undef DLDPQ_KEY_T
#undef DLDPQ_VAL_T


#define DLRAND_PREFIX sint
#define DLRAND_TYPE_T sint_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX sint
#define DLSTATS_TYPE_T sint_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLLIST_PREFIX sint
#define DLLIST_TYPE_T sint_t
#define DLLIST_LINKED
#include "dllist_headers.h"
#undef DLLIST_LINKED
#undef DLLIST_PREFIX
#undef DLLIST_TYPE_T


#define DLSORT_PREFIX sint
#define DLSORT_TYPE_T sint_t
#include "dlsort_headers.h"
#undef DLSORT_PREFIX
#undef DLSORT_TYPE_T


#define DLHEAP_PREFIX sint
#define DLHEAP_TYPE_T sint_t
#include "dlheap_headers.h"
#undef DLHEAP_PREFIX
#undef DLHEAP_TYPE_T


#define DLTREE_PREFIX sint
#define DLTREE_KEY_T char*
#define DLTREE_VAL_T sint_t
#include "dltree_headers.h"
#undef DLTREE_KEY_T
#undef DLTREE_VAL_T
#undef DLTREE_PREFIX




/* uint_t */
#define DLMEM_PREFIX uint 
#define DLMEM_TYPE_T uint_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLMATH_PREFIX uint
#define DLMATH_TYPE_T uint_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


DL_MK_OMP_REDUCTION_HEADERS(sint,sint_t)
DL_MK_SORTKV_HEADERS(sint,sint_t,sint_t)



#define TESTEQUALS(a,b,fmt) \
  do { \
    if ((a) != (b)) { \
      eprintf("[TEST] " #a " (" fmt ") != " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define OMPTESTEQUALS(a,b,fmt,r) \
  do { \
    if ((a) != (b)) { \
      eprintf("[TEST]<tid=%d> " #a " (" fmt ") != " #b " (" fmt \
          ") at %s:%d\n",omp_get_thread_num(),a,b,__FILE__,__LINE__); \
      _Pragma("omp atomic") \
      ++(r); \
    } \
  } while(0)

#define TESTLESSTHANOREQUAL(a,b,fmt) \
  do { \
    if ((a) > (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! <= " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTGREATERTHANOREQUAL(a,b,fmt) \
  do { \
    if ((a) < (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! >= " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTLESSTHAN(a,b,fmt) \
  do { \
    if ((a) >= (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! < " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTGREATERTHAN(a,b,fmt) \
  do { \
    if ((a) <= (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! > " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTTRUE(a) \
  do { \
    if (!(a)) { \
      eprintf("[TEST] " #a " evaluated to false at %s:%d\n",__FILE__, \
          __LINE__); \
      return 1; \
    } \
  } while (0)

#define OMPTESTTRUE(a,r) \
  do { \
    if (!(a)) { \
      eprintf("[TEST] " #a " evaluated to false at %s:%d\n",__FILE__, \
          __LINE__); \
      _Pragma("omp atomic") \
      ++(r); \
    } \
  } while (0)



sint_t test(void);

#endif
