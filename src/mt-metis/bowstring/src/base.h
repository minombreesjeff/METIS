/**
 * @file base.h
 * @brief Library external header file
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Dominique LaSalle
 * @version 1
 * @date 2013-04-09
 */




#ifndef BOWSTRING_BASE_H
#define BOWSTRING_BASE_H




#include <domlib.h>
#include <stdlib.h> 
#include <stdio.h>
#include <stdint.h>
#include <bowstring.h>
#include "strings.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/**
 * @brief Flags specifying what weights are in a file
 */
typedef enum bowstring_wgt_flag_t {
  ADJWGT_FLAG = 1,
  VWGT_FLAG = 1 << 1,
  IADJWGT_FLAG = 1 << 2
} bowstring_wgt_flag_t;


/* rename bowstring types to short names */
#define vtx_t bowstring_vtx_t
#define adj_t bowstring_adj_t
#define wgt_t bowstring_wgt_t
#define vlbl_t bowstring_vlbl_t
#define elbl_t bowstring_elbl_t
#define coord_t bowstring_coord_t


/* rename dynamic graph types to short names */
#define dynnode_t bowstring_dynnode_t
#define dyngraph_t bowstring_dyngraph_t




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


/* wgt_t */
#define DLMEM_PREFIX wgt
#define DLMEM_TYPE_T wgt_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX wgt
#define DLMATH_TYPE_T wgt_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX wgt
#define DLRAND_TYPE_T wgt_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


/* vtx_t */
#define DLMEM_PREFIX vtx
#define DLMEM_TYPE_T vtx_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLRAND_PREFIX vtx
#define DLRAND_TYPE_T vtx_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLMATH_PREFIX vtx
#define DLMATH_TYPE_T vtx_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* adj_t */
#define DLMEM_PREFIX adj
#define DLMEM_TYPE_T adj_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLRAND_PREFIX adj
#define DLRAND_TYPE_T adj_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLMATH_PREFIX adj
#define DLMATH_TYPE_T adj_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* vlbl_t */
#define DLMEM_PREFIX vlbl
#define DLMEM_TYPE_T vlbl_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX vlbl
#define DLMATH_TYPE_T vlbl_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* elbl_t */
#define DLMEM_PREFIX elbl
#define DLMEM_TYPE_T elbl_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


/* coord_t */
#define DLMEM_PREFIX coord
#define DLMEM_TYPE_T coord_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLRAND_PREFIX coord
#define DLRAND_TYPE_T coord_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const vtx_t NULL_VTX = (vtx_t)-1;
static const adj_t NULL_ADJ = (adj_t)-1;




#endif
