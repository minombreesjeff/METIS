#include "test.h"

/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#include "dlmem_funcs.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX real
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLPQ_PREFIX rs
#define DLPQ_KEY_T real_t
#define DLPQ_VAL_T sint_t
#include "dlpq_funcs.h"
#undef DLPQ_PREFIX
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T




/* sint_t */
#define DLMEM_PREFIX sint 
#define DLMEM_TYPE_T sint_t
#include "dlmem_funcs.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX sint
#define DLMATH_TYPE_T sint_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLBUFFER_PREFIX sint
#define DLBUFFER_TYPE_T sint_t
#include "dlbuffer_funcs.h"
#undef DLBUFFER_PREFIX
#undef DLBUFFER_TYPE_T


#define DLHT_PREFIX ss
#define DLHT_KEY_T sint_t
#define DLHT_VAL_T sint_t
#include "dlht_funcs.h"
#undef DLHT_PREFIX
#undef DLHT_KEY_T
#undef DLHT_VAL_T


#define DLCB_PREFIX sint
#define DLCB_TYPE_T sint_t
#include "dlcb_funcs.h"
#undef DLCB_PREFIX
#undef DLCB_TYPE_T


#define DLISET_PREFIX sint
#define DLISET_TYPE_T sint_t
#include "dliset_funcs.h"
#undef DLISET_PREFIX
#undef DLISET_TYPE_T


#define DLDJSET_PREFIX sint
#define DLDJSET_TYPE_T sint_t
#include "dldjset_funcs.h"
#undef DLDJSET_PREFIX
#undef DLDJSET_TYPE_T


#define DLDPQ_PREFIX sint
#define DLDPQ_KEY_T sint_t
#define DLDPQ_VAL_T sint_t
#include "dldpq_funcs.h"
#undef DLDPQ_PREFIX
#undef DLDPQ_KEY_T
#undef DLDPQ_VAL_T


#define DLSTATS_PREFIX sint 
#define DLSTATS_TYPE_T sint_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLSORT_PREFIX sint
#define DLSORT_TYPE_T sint_t
#define DLSORT_DLSIGN DLSIGN_INTEGRAL
#include "dlsort_funcs.h"
#undef DLSORT_DLSIGN
#undef DLSORT_PREFIX
#undef DLSORT_TYPE_T


#define DLLIST_PREFIX sint
#define DLLIST_TYPE_T sint_t
#define DLLIST_LINKED
#include "dllist_funcs.h"
#undef DLLIST_LINKED
#undef DLLIST_PREFIX
#undef DLLIST_TYPE_T


#define DLHEAP_PREFIX sint
#define DLHEAP_TYPE_T sint_t
#define DLHEAP_MAX
#include "dlheap_funcs.h"
#undef DLHEAP_MAX
#undef DLHEAP_PREFIX
#undef DLHEAP_TYPE_T


#define DLRAND_PREFIX sint
#define DLRAND_TYPE_T sint_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLTREE_PREFIX sint
#define DLTREE_KEY_T char*
#define DLTREE_VAL_T sint_t
#include "dltree_funcs.h"
#undef DLTREE_KEY_T
#undef DLTREE_VAL_T
#undef DLTREE_PREFIX




/* uint_t */
#define DLMEM_PREFIX uint 
#define DLMEM_TYPE_T uint_t
#include "dlmem_funcs.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX uint
#define DLMATH_TYPE_T uint_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


DL_MK_OMP_REDUCTION_FUNCS(sint,sint_t)
DL_MK_SORTKV_FUNCS(sint,sint_t,sint_t)


sint_t main(const int argc, const char ** argv) 
{
  dl_timer_t tmr;
  dl_init_timer(&tmr);
  dl_start_timer(&tmr);
  int rv = test();
  double time = dl_poll_timer(&tmr);

  if (rv == 0) {
    printf("PASSED %7.3f s\n",time);
  } else { 
    printf("FAILED %7.3f s\n",time);
  }

  return rv; 
}
