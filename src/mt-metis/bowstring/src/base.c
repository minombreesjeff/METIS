/**
 * @file base.c
 * @brief Base library functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Dominique LaSalle
 * @version 1
 * @date 2013-04-09
 */




#ifndef BOWSTRING_BASE_C
#define BOWSTRING_BASE_C




#include "base.h"




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#ifndef BOWSTRING_TYPES_DEFINED


/* wgt_t */
#define DLMEM_PREFIX wgt
#define DLMEM_TYPE_T wgt_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLMATH_PREFIX wgt 
#define DLMATH_TYPE_T wgt_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T

#define DLRAND_PREFIX wgt
#define DLRAND_TYPE_T wgt_t
#define DLRAND_DLTYPE DLTYPE_FLOAT
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


/* vtx_t */
#define DLMEM_PREFIX vtx
#define DLMEM_TYPE_T vtx_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLRAND_PREFIX vtx
#define DLRAND_TYPE_T vtx_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T

#define DLMATH_PREFIX vtx
#define DLMATH_TYPE_T vtx_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* adj_t */
#define DLMEM_PREFIX adj
#define DLMEM_TYPE_T adj_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLRAND_PREFIX adj
#define DLRAND_TYPE_T adj_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T

#define DLMATH_PREFIX adj
#define DLMATH_TYPE_T adj_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#endif




#ifndef BOWSTRING_LABELS_DEFINED


/* vlbl_t */
#define DLMEM_PREFIX vlbl
#define DLMEM_TYPE_T vlbl_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLMATH_PREFIX vlbl
#define DLMATH_TYPE_T vlbl_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* elbl_t */
#define DLMEM_PREFIX elbl
#define DLMEM_TYPE_T elbl_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#endif



#ifndef BOWSTRING_COORDS_DEFINED


/* coord_t */
#define DLMEM_PREFIX coord
#define DLMEM_TYPE_T coord_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLRAND_PREFIX coord
#define DLRAND_TYPE_T coord_t
#define DLRAND_DLTYPE DLTYPE_FLOAT
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#endif


#endif
