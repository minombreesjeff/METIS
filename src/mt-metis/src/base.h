/**
 * @file base.h
 * @brief Base types etc.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */




#ifndef MTMETIS_BASE_H
#define MTMETIS_BASE_H




/******************************************************************************
* EXTERNAL INCLUDES ***********************************************************
******************************************************************************/


#define _POSIX_SOURCE 1

#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <mtmetis.h>
#include <domlib.h>




/******************************************************************************
* MACRO INCLUDES **************************************************************
******************************************************************************/


#include "strings.h"
#include "macros.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


#ifdef MTMETIS_64BIT_THREADS
typedef uint64_t mtmetis_tid_t;
#else
typedef uint32_t mtmetis_tid_t;
#endif


#ifdef MTMETIS_DOUBLE_REAL
typedef double mtmetis_real_t;
#else
typedef float mtmetis_real_t;
#endif

typedef int64_t mtmetis_twgt_t;
typedef uint16_t mtmetis_offset_t;


/* rename mtmetis.h types for internal use */
#define vtx_t mtmetis_vtx_t
#define adj_t mtmetis_adj_t
#define wgt_t mtmetis_wgt_t
#define twgt_t mtmetis_twgt_t
#define pid_t mtmetis_pid_t
#define tid_t mtmetis_tid_t
#define real_t mtmetis_real_t
#define offset_t mtmetis_offset_t




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


/* macros */
#define DEF_NULL_VTX ((vtx_t)-1)
#define DEF_NULL_ADJ ((adj_t)-1)
#define DEF_NULL_PID ((pid_t)-1)
#define DEF_NULL_TID ((tid_t)-1)
#define DEF_NULL_WGT ((wgt_t)-1)
#define DEF_NULL_OFFSET ((offset_t)-1)


/* type null values */
static const vtx_t NULL_VTX = DEF_NULL_VTX;
static const wgt_t NULL_WGT = DEF_NULL_WGT;
static const adj_t NULL_ADJ = DEF_NULL_ADJ;
static const pid_t NULL_PID = DEF_NULL_PID;
static const tid_t NULL_TID = DEF_NULL_TID;
static const offset_t NULL_OFFSET = DEF_NULL_OFFSET;


/* thread specific constants */
static const vtx_t BLOCKSIZE = 0x1000;
static const int BLOCKSHIFT = 12;
static const vtx_t BLOCKMASK = 0x0FFF;




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


/* vtx_t */
#define DLMEM_PREFIX vtx
#define DLMEM_TYPE_T vtx_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX vtx
#define DLMATH_TYPE_T vtx_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX vtx
#define DLRAND_TYPE_T vtx_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#define DLRAND_STATIC
#include "dlrand_headers.h"
#undef DLRAND_STATIC
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX vtx
#define DLSTATS_TYPE_T vtx_t
#define DLSTATS_STATIC
#include "dlstats_headers.h"
#undef DLSTATS_STATIC
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLISET_PREFIX vtx
#define DLISET_TYPE_T vtx_t
#define DLISET_STATIC
#include "dliset_headers.h"
#undef DLISET_STATIC
#undef DLISET_TYPE_T
#undef DLISET_PREFIX


#define DLTHREAD_PREFIX vtx
#define DLTHREAD_TYPE_T vtx_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX


#define DLSORT_PREFIX vtx
#define DLSORT_TYPE_T vtx_t
#define DLSORT_STATIC
#include "dlsort_headers.h"
#undef DLSORT_STATIC
#undef DLSORT_TYPE_T
#undef DLSORT_PREFIX




/* adj_t */
#define DLMEM_PREFIX adj
#define DLMEM_TYPE_T adj_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX adj
#define DLMATH_TYPE_T adj_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX adj
#define DLRAND_TYPE_T adj_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#define DLRAND_STATIC
#include "dlrand_headers.h"
#undef DLRAND_STATIC
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX adj
#define DLSTATS_TYPE_T adj_t
#define DLSTATS_STATIC
#include "dlstats_headers.h"
#undef DLSTATS_STATIC
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLTHREAD_PREFIX adj
#define DLTHREAD_TYPE_T adj_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX




/* pid_t */
#define DLMEM_PREFIX pid
#define DLMEM_TYPE_T pid_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX pid
#define DLMATH_TYPE_T pid_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX pid
#define DLRAND_TYPE_T pid_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#define DLRAND_STATIC
#include "dlrand_headers.h"
#undef DLRAND_STATIC
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX pid
#define DLSTATS_TYPE_T pid_t
#define DLSTATS_STATIC
#include "dlstats_headers.h"
#undef DLSTATS_STATIC
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


/* tid_t */
#define DLMEM_PREFIX tid
#define DLMEM_TYPE_T tid_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX tid
#define DLMATH_TYPE_T tid_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* wgt_t */
#define DLMEM_PREFIX wgt
#define DLMEM_TYPE_T wgt_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX wgt
#define DLMATH_TYPE_T wgt_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX wgt
#define DLRAND_TYPE_T wgt_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#define DLRAND_STATIC
#include "dlrand_headers.h"
#undef DLRAND_STATIC
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX wgt
#define DLSTATS_TYPE_T wgt_t
#define DLSTATS_STATIC
#include "dlstats_headers.h"
#undef DLSTATS_STATIC
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLTHREAD_PREFIX wgt
#define DLTHREAD_TYPE_T wgt_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX


/* twgt_t */
#define DLTHREAD_PREFIX twgt
#define DLTHREAD_TYPE_T twgt_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX


#define DLMATH_PREFIX twgt
#define DLMATH_TYPE_T twgt_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX real
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLTHREAD_PREFIX double
#define DLTHREAD_TYPE_T double
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX


/* int */

#define DLTHREAD_PREFIX int
#define DLTHREAD_TYPE_T int
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX


/* offset_t */

#define DLMEM_PREFIX offset
#define DLMEM_TYPE_T offset_t
#define DLMEM_STATIC 1
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#ifdef DEBUG
  #define par_dprintf(...) \
    do { \
      _Pragma("omp master") \
      { \
        dprintf( __VA_ARGS__ ); \
      } \
    } while(0)
#else
  #define par_dprintf(...)
#endif

#define par_vprintf(...) \
  do { \
    _Pragma("omp master") \
    { \
      vprintf( __VA_ARGS__ ); \
    } \
  } while(0)




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


static inline vtx_t gvtx_to_lvtx(
    vtx_t const v, 
    vtx_t const mask)
{
  DL_ASSERT(mask > 0,"The mask is set to 0!\n");
  DL_ASSERT(v > mask,"Global vertex number is smaller than mask (gvtx = %"
      PF_VTX_T", mask = %"PF_VTX_T")\n",v,mask);
  return v & mask;
}


static inline vtx_t lvtx_to_gvtx(
    vtx_t const v, 
    tid_t const t, 
    int const shift)
{
  DL_ASSERT(shift > 0,"The mask size is set to 0!\n");
  DL_ASSERT(v < (vtx_t)(1 << shift),"Local vertex number is greater than "
      "shift (lvtx = %"PF_VTX_T", shift = %d)\n",v,shift);
  return ((t+1) << shift) | v;
}


static inline tid_t gvtx_to_tid(
    vtx_t const v, 
    int const shift)
{
  DL_ASSERT(shift > 0,"The shift size is set to %d!\n",shift);
  DL_ASSERT(v >= (vtx_t)(1 << shift),"Global vertex number is too small "
      "(gvtx = %"PF_VTX_T", shift = %d)\n",v,shift);
  return (v >> shift)-1;
}


static inline vtx_t max_gvtx(
    int const shift, 
    tid_t const nthreads) 
{
  return (vtx_t)(1 << shift)*(nthreads+1);
}


static inline int is_bnd(
    wgt_t const id,
    wgt_t const ed,
    int const greedy)
{
  if (greedy) {
    return ed >= id;
  } else {
    return (ed > 0) || (id == 0);
  }
}


/* avoid having to pass each element */
#define gvtx_to_lvtx(v,dist) gvtx_to_lvtx(v,(dist).mask)
#define lvtx_to_gvtx(v,t,dist) lvtx_to_gvtx(v,t,(dist).shift)
#define gvtx_to_tid(v,dist) gvtx_to_tid(v,(dist).shift)
#define max_gvtx(graph) max_gvtx((graph)->dist.shift,(graph)->dist.nthreads)






#endif
